"""
Finding the upstream outputs.

Linked layout: each input port arrives in INPUT_DIR/<port>/ (placement: subfolder), its entries
often symlinks to the upstream outputs. Flat layout (no links, everything side by side): files are
found by name anywhere under INPUT_DIR. Files come from the shared walker (processor/safe_inputs.py:
links that leave the input are skipped) and are opened only through it, capped; a file "next to"
another one is used only when the walker accepted it. An upstream folder is recognised by the files
only its app writes:

  bulk-de       a folder with design.json and summary.json (one per LABEL; several may be linked)
  rnaseq-quant  a samples.csv with reads_in and mapping_rate columns (+ qc/<sample>.fastp.json)

Every provenance.json with an "app" is read for the methods and citations. The per-gene tables are
read only for the top genes; nothing per read is ever opened.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from processor.safe_inputs import InputDir, UnsafeInput

log = logging.getLogger("bulk-report")
MAX_JSON = 20 << 20
MAX_CSV = 200 << 20
MAX_PNG = 10 << 20


class Files:
    """The accepted files of INPUT_DIR and capped readers that return None (logged) on any problem."""

    def __init__(self, input_dir: Path):
        self.inp = InputDir(input_dir)
        self.all = self.inp.walk()
        self._set = {str(p) for p in self.all}

    def rel(self, path: Path) -> str:
        return self.inp.rel(path)

    def get(self, folder: Path, name: str) -> Path | None:
        p = folder / name
        return p if str(p) in self._set else None

    def json(self, path: Path | None) -> dict[str, Any] | None:
        if path is None:
            return None
        try:
            value = self.inp.read_json(path, MAX_JSON)
        except (UnsafeInput, OSError, ValueError) as exc:
            log.warning("cannot read %s: %s", self.rel(path), exc)
            return None
        return value if isinstance(value, dict) else None

    def csv(self, path: Path | None, *cols: str, **kwargs) -> pd.DataFrame | None:
        if path is None:
            return None
        try:
            df = self.inp.read_csv(path, MAX_CSV, **kwargs)
        except (UnsafeInput, OSError, ValueError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
            log.warning("cannot read %s: %s", self.rel(path), exc)
            return None
        if any(c not in df.columns for c in cols):
            log.warning("%s: not the expected table; ignored", self.rel(path))
            return None
        return df

    def png(self, path: Path | None) -> str | None:
        return self.inp.inline_png(path, MAX_PNG) if path is not None else None


@dataclass
class DeRun:
    dir: Path
    label: str
    summary: dict[str, Any]
    design: dict[str, Any]
    de: dict[str, Any]
    samples: pd.DataFrame | None
    pca: pd.DataFrame | None
    files: Files

    def file(self, name: str) -> Path | None:
        return self.files.get(self.dir, name)


@dataclass
class Quant:
    dir: Path
    samples: pd.DataFrame
    fastp: dict[str, dict[str, Any]] = field(default_factory=dict)     # sample_id -> fastp summary


@dataclass
class Inputs:
    files: Files
    de_runs: list[DeRun] = field(default_factory=list)
    quant: Quant | None = None
    provenances: list[tuple[Path, dict[str, Any]]] = field(default_factory=list)

    def rel(self, path: Path) -> str:
        return self.files.rel(path)


def _port_files(files: Files, port: str) -> list[Path]:
    inside = [f for f in files.all if files.rel(f).startswith(f"{port}/")]
    return inside or files.all


def gather(input_dir: Path) -> Inputs:
    files = Files(input_dir)
    inp = Inputs(files)
    seen_labels: set[str] = set()
    for f in _port_files(files, "de"):
        if f.name != "design.json":
            continue
        summary, design = files.json(files.get(f.parent, "summary.json")), files.json(f)
        if not summary or not design or not isinstance(summary.get("contrasts"), list):
            continue
        label = str(summary.get("label") or f.parent.name)
        if label in seen_labels:
            label = f"{label} ({files.rel(f.parent)})"
        seen_labels.add(label)
        inp.de_runs.append(DeRun(f.parent, label, summary, design, files.json(files.get(f.parent, "de.json")) or {},
                                 files.csv(files.get(f.parent, "samples.csv"), "sample_id", dtype={"sample_id": str}),
                                 files.csv(files.get(f.parent, "pca.csv"), "sample_id", "PC1", dtype={"sample_id": str}),
                                 files))
    for f in _port_files(files, "counts"):
        if f.name != "samples.csv":
            continue
        df = files.csv(f, "sample_id", "reads_in", "mapping_rate", dtype={"sample_id": str})
        if df is not None:
            q = Quant(f.parent, df)
            for j in files.all:
                if j.parent == f.parent / "qc" and j.name.endswith(".fastp.json"):
                    data = files.json(j)
                    if data and isinstance(data.get("summary"), dict):
                        q.fastp[j.name.removesuffix(".fastp.json")] = data["summary"]
            inp.quant = q
            break
    for f in files.all:
        if f.name == "provenance.json":
            prov = files.json(f)
            if prov and prov.get("app"):
                inp.provenances.append((f, prov))
    return inp
