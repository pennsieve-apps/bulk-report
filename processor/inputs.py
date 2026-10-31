"""
Finding the upstream outputs.

Linked layout: each input port arrives in INPUT_DIR/<port>/ (placement: subfolder).
Flat layout (no links, everything side by side): files are found by name anywhere under
INPUT_DIR. Either way an upstream folder is recognised by the files only its app writes:

  bulk-de       a folder with design.json and summary.json (one per LABEL; several may be linked)
  rnaseq-quant  a samples.csv with reads_in and mapping_rate columns (+ qc/<sample>.fastp.json)

Every provenance.json with an "app" is read for the methods and citations. The walk
follows symlinks (INPUT_DIR entries are symlinks) but not links below an entry that leave
its tree, visits a folder once, and skips files over the size caps. The per-gene tables
are read only for the top genes; nothing per-read is ever opened.
"""
from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

log = logging.getLogger("bulk-report")
MAX_JSON = 20 << 20
MAX_CSV = 512 << 20


def _inside(path: str, root: str) -> bool:
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


def walk(root: Path) -> list[Path]:
    found: list[Path] = []
    root_real = os.path.realpath(root)
    for entry in sorted(root.iterdir()) if root.is_dir() else []:
        if entry.name.startswith("."):
            continue
        base = os.path.realpath(entry)
        if entry.is_file():
            found.append(entry)
            continue
        seen: set[str] = set()
        for dirpath, dirnames, filenames in os.walk(entry, followlinks=True):
            real = os.path.realpath(dirpath)
            if real in seen or not (_inside(real, base) or _inside(real, root_real)):
                dirnames[:] = []
                continue
            seen.add(real)
            dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
            for name in sorted(filenames):
                p = Path(dirpath) / name
                rp = os.path.realpath(p)
                if not name.startswith(".") and os.path.isfile(rp) and (_inside(rp, base) or _inside(rp, root_real)):
                    found.append(p)
    return found


def read_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not path.is_file() or path.stat().st_size > MAX_JSON:
        return None
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else None
    except (OSError, ValueError) as exc:
        log.warning("cannot read %s: %s", path, exc)
        return None


def read_csv(path: Path | None, *cols: str, **kwargs) -> pd.DataFrame | None:
    if path is None or not path.is_file() or path.stat().st_size > MAX_CSV:
        return None
    try:
        df = pd.read_csv(path, **kwargs)
    except (OSError, ValueError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        log.warning("cannot read %s: %s", path, exc)
        return None
    if any(c not in df.columns for c in cols):
        log.warning("%s: not the expected table; ignored", path)
        return None
    return df


@dataclass
class DeRun:
    dir: Path
    label: str
    summary: dict[str, Any]
    design: dict[str, Any]
    de: dict[str, Any]
    samples: pd.DataFrame | None
    pca: pd.DataFrame | None

    def file(self, name: str) -> Path | None:
        p = self.dir / name
        return p if p.is_file() else None


@dataclass
class Quant:
    dir: Path
    samples: pd.DataFrame
    fastp: dict[str, dict[str, Any]] = field(default_factory=dict)     # sample_id -> fastp summary


@dataclass
class Inputs:
    input_dir: Path
    de_runs: list[DeRun] = field(default_factory=list)
    quant: Quant | None = None
    provenances: list[tuple[Path, dict[str, Any]]] = field(default_factory=list)

    def rel(self, path: Path) -> str:
        try:
            return str(Path(path).relative_to(self.input_dir))
        except ValueError:
            return str(path)


def _port_files(files: list[Path], input_dir: Path, port: str) -> list[Path]:
    d = input_dir / port
    inside = [f for f in files if _inside(str(f), str(d))] if d.is_dir() else []
    return inside or files


def gather(input_dir: Path) -> Inputs:
    files = walk(input_dir)
    inp = Inputs(input_dir)
    seen_labels: set[str] = set()
    for f in _port_files(files, input_dir, "de"):
        if f.name != "design.json" or not (f.parent / "summary.json").is_file():
            continue
        summary, design = read_json(f.parent / "summary.json"), read_json(f)
        if not summary or not design or "contrasts" not in summary:
            continue
        label = str(summary.get("label") or f.parent.name)
        if label in seen_labels:
            label = f"{label} ({inp.rel(f.parent)})"
        seen_labels.add(label)
        inp.de_runs.append(DeRun(f.parent, label, summary, design, read_json(f.parent / "de.json") or {},
                                 read_csv(f.parent / "samples.csv", "sample_id", dtype={"sample_id": str}),
                                 read_csv(f.parent / "pca.csv", "sample_id", "PC1", dtype={"sample_id": str})))
    for f in _port_files(files, input_dir, "counts"):
        if f.name != "samples.csv":
            continue
        df = read_csv(f, "sample_id", "reads_in", "mapping_rate", dtype={"sample_id": str})
        if df is not None:
            q = Quant(f.parent, df)
            for j in sorted((f.parent / "qc").glob("*.fastp.json")) if (f.parent / "qc").is_dir() else []:
                data = read_json(j)
                if data and isinstance(data.get("summary"), dict):
                    q.fastp[j.name.removesuffix(".fastp.json")] = data["summary"]
            inp.quant = q
            break
    for f in files:
        if f.name == "provenance.json":
            prov = read_json(f)
            if prov and prov.get("app"):
                inp.provenances.append((f, prov))
    return inp
