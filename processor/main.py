"""
bulk-report: one offline HTML page from the bulk RNA-seq workflow's outputs.

Reads bulk-de's <LABEL>/ folders (port `de`, one or more) and, when linked,
rnaseq-quant's counts output (port `counts`: samples.csv with read and mapping
numbers, qc/*.fastp.json), and writes to OUTPUT_DIR:

  index.html       design, sample QC, PCA, per contrast volcano / MA / top genes / GSEA,
                   methods (every step's version and parameters), citations gathered from
                   every provenance.json, data credits. Figures are embedded (data URIs): the
                   page is one file, with no scripts and no external resources.
  summary.json     the page's headline numbers, observations and caveats
  citations.csv    the numbered citations
  provenance.json
"""
from __future__ import annotations

import base64
import csv
import json
import logging
import os
import platform
import sys
import tempfile
import time
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("bulk-report")
log.setLevel(logging.INFO)
for noisy in ("matplotlib", "PIL", "fontTools"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "bulk-report"
REPO_ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = Path(__file__).parent / "templates"

DEFAULTS = {"REPORT_TITLE": "Bulk RNA-seq Report", "FDR": "0.05", "MAX_GENES": "15", "DATA_CREDIT": ""}
PARAMETERS = tuple(DEFAULTS)


class ParameterError(ValueError):
    pass


def read_parameters(env=os.environ) -> dict:
    raw = {k: (str(env.get(k)).strip() if env.get(k) not in (None, "") else v) for k, v in DEFAULTS.items()}
    try:
        p = {"REPORT_TITLE": raw["REPORT_TITLE"][:200], "FDR": float(raw["FDR"]), "MAX_GENES": int(raw["MAX_GENES"]),
             "DATA_CREDIT": raw["DATA_CREDIT"][:2000]}
    except ValueError as exc:
        raise ParameterError(f"bad numeric parameter: {exc}") from exc
    if not 0 < p["FDR"] < 1:
        raise ParameterError("FDR must be in (0, 1)")
    if not 1 <= p["MAX_GENES"] <= 200:
        raise ParameterError("MAX_GENES must be in 1..200")
    return p


def app_version() -> str:
    try:
        return (REPO_ROOT / "VERSION").read_text().strip()
    except OSError:
        return "unknown"


def tool_versions() -> dict:
    out = {"python": platform.python_version()}
    for pkg in ("pandas", "matplotlib", "jinja2"):
        try:
            out[pkg] = version(pkg)
        except PackageNotFoundError:
            out[pkg] = None
    return out


def qc_figure(F, qc) -> str:
    """The report's own QC figure as a data URI (drawn into a temporary folder of its own)."""
    with tempfile.TemporaryDirectory() as tmp:
        return "data:image/png;base64," + base64.b64encode(F.qc_bars(qc, Path(tmp) / "qc.png").read_bytes()).decode()


def own_citations() -> list[dict]:
    from processor.sections import _references, _ref_text

    refs = _references()
    return [{"doi": refs[k]["doi"], "text": _ref_text(refs[k])} | ({"pmid": refs[k]["pmid"]} if refs[k].get("pmid") else {})
            for k in ("matplotlib",)]


def peak_rss_mb():
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return round(peak / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)
    except (ImportError, OSError):
        return None


def main() -> int:
    started = datetime.now(timezone.utc)
    t0 = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 2
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    try:
        params = read_parameters()
    except ParameterError as exc:
        log.error("%s", exc)
        return 2

    from jinja2 import Environment, FileSystemLoader, select_autoescape

    from processor import figures as F
    from processor import sections as S
    from processor.inputs import gather

    inp = gather(input_dir)
    log.info("found %d bulk-de run(s) (%s), rnaseq-quant QC: %s, %d provenance file(s)", len(inp.de_runs),
             ", ".join(r.label for r in inp.de_runs) or "none", "yes" if inp.quant else "no", len(inp.provenances))
    if not inp.de_runs and inp.quant is None:
        log.error("no bulk-de output (design.json + summary.json) or rnaseq-quant samples.csv under %s", input_dir)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    errors: list[dict] = []

    def section(name: str, fn, default=None):
        """Each section on its own: an unexpected upstream file costs that section, not the page."""
        try:
            return fn()
        except Exception as exc:
            log.exception("%s: %s", name, exc)
            errors.append({"section": name, "error": f"{type(exc).__name__}: {exc}"})
            return default

    qc = section("sample QC", lambda: S.sample_qc(inp))
    qc_png = section("sample QC figure", lambda: qc_figure(F, qc)) if qc else None
    runs = [x for x in (section(f"differential expression {r.label}",
                                lambda r=r: S.de_run(r, params["FDR"], params["MAX_GENES"])) for r in inp.de_runs) if x]
    own = {"app": APP_NAME, "version": app_version(), "tools": tool_versions(),
           "parameters": {k: v for k, v in params.items() if k != "DATA_CREDIT"}, "citations": own_citations()}
    meth = section("methods and citations", lambda: S.methods(inp, own), {"steps": [], "references": []})
    credit = section("data credits", lambda: S.credits(inp, params["DATA_CREDIT"]), [])
    generated = started.strftime("%Y-%m-%d %H:%M UTC")

    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=select_autoescape(["html", "j2"]),
                      trim_blocks=True, lstrip_blocks=True)
    env.filters.update(num=S.num, pval=S.pval, pct=S.pct, gsea_p=S.gsea_p)
    html = env.get_template("index.html.j2").render(title=params["REPORT_TITLE"], generated_at=generated, qc=qc,
                                                    qc_png=qc_png, runs=runs, methods=meth, credit=credit,
                                                    fdr=params["FDR"], errors=errors)
    (output_dir / "index.html").write_text(html, encoding="utf-8")
    with open(output_dir / "citations.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["number", "citation", "doi", "pmid", "verified"])
        for i, r in enumerate(meth["references"], 1):
            w.writerow([i, r["text"], r["doi"], r["pmid"], "yes" if r.get("verified") else "no"])
    summary = {
        "app": APP_NAME, "version": own["version"], "generated": generated, "title": params["REPORT_TITLE"],
        "fdr": params["FDR"],
        "sample_qc": {k: qc[k] for k in ("source", "observations", "caveats", "n_flagged")} if qc else None,
        "differential_expression": [{
            "label": r["label"], "formula": r["formula"], "n_samples": len(r["samples"]),
            "observations": r["observations"], "caveats": r["caveats"],
            "contrasts": [{k: c.get(k) for k in ("contrast", "n_up", "n_down", "n_significant")}
                          | {"top_up": [g["gene_name"] or g["gene_id"] for g in c["up"][:10]],
                             "top_down": [g["gene_name"] or g["gene_id"] for g in c["down"][:10]],
                             "gsea": [{k: g[k] for k in ("collection", "n_sets", "n_up", "n_down")} for g in c["gsea"]]}
                          for c in r["contrasts"]]} for r in runs],
        "steps": [{k: s[k] for k in ("app", "version", "where", "citations")} for s in meth["steps"]],
        "citations": [{"number": i, **r} for i, r in enumerate(meth["references"], 1)],
        "data": credit, "errors": errors,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    prov = {**own, "inputs": {"de": [inp.rel(r.dir) for r in inp.de_runs],
                              "counts": inp.rel(inp.quant.dir) if inp.quant else None,
                              "provenance": [inp.rel(p) for p, _ in inp.provenances]},
            "started": started.isoformat(), "finished": datetime.now(timezone.utc).isoformat(),
            "samples": [{"sample_id": r["sample_id"], "flags": r["flags"]} for r in (qc or {}).get("rows", [])],
            "peak_rss_mb": peak_rss_mb(), "timings_seconds": {"total": round(time.monotonic() - t0, 2)},
            "outputs": ["index.html", "summary.json", "citations.csv", "provenance.json"], "errors": errors}
    (output_dir / "provenance.json").write_text(json.dumps(prov, indent=2, default=str) + "\n")
    log.info("done: %d run(s), %d contrast(s), %d citations, %.1fs -> %s", len(runs),
             sum(len(r["contrasts"]) for r in runs), len(meth["references"]), time.monotonic() - t0,
             output_dir / "index.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
