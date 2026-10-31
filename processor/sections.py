"""
The report's sections, from the upstream outputs: plain dicts for the template, with
observations (numbers the data shows) and caveats (what the method cannot tell).
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from processor.inputs import DeRun, Inputs

LOW_MAPPING = 0.5           # salmon mapping rate below which a library is flagged
LOW_READS = 5_000_000       # reads after trimming below which a library is flagged
OUTLIER_FOLD = 3.0          # a library size this many times below the median is flagged

REFERENCES = Path(__file__).parent / "references.json"


def num(v, d: int = 0) -> str:
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "–"
    return f"{float(v):,.{d}f}"


def pval(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "–"
    v = float(v)
    return "< 1e-300" if v == 0 else (f"{v:.2g}" if v < 0.001 else f"{v:.3f}")


def gsea_p(v, permutations: int = 1000) -> str:
    """GSEA p/FDR from permutations: 0 means smaller than 1/(permutations + 1)."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "–"
    return f"< 1/{permutations + 1}" if float(v) == 0 else pval(v)


def pct(v, d: int = 1) -> str:
    return "–" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{float(v) * 100:.{d}f}%"


# ---- sample QC ---------------------------------------------------------------------------------

def sample_qc(inp: Inputs) -> dict[str, Any] | None:
    q = inp.quant
    if q is None:
        run = next((r for r in inp.de_runs if r.samples is not None), None)
        if run is None:
            return None
        s = run.samples
        med = float(s["library_size"].median())
        rows = []
        for _, r in s.iterrows():
            flags = []
            if r["library_size"] * OUTLIER_FOLD < med:
                flags.append(f"library {OUTLIER_FOLD:g}x below the median")
            rows.append({"sample_id": r["sample_id"], "library_size": int(r["library_size"]),
                         "n_genes_detected": int(r["n_genes_detected"]), "size_factor": r.get("size_factor"),
                         "used": bool(r["used"]), "flags": flags})
        obs = [f"{len(rows)} libraries from {run.de.get('layout', 'a count table')}; counts per library "
               f"{num(s['library_size'].min())} to {num(s['library_size'].max())} (median {num(med)}); "
               f"{num(s['n_genes_detected'].median())} genes detected (median)."]
        return {"source": "counts", "rows": rows, "observations": obs,
                "caveats": ["No read-level QC (trimming, mapping rate) was linked: the counts were deposited or "
                            "computed elsewhere. Library size and genes detected are the only checks here."],
                "n_flagged": sum(1 for r in rows if r["flags"])}
    s = q.samples.copy()
    for c in ("reads_in", "reads_after_trim", "mapping_rate"):
        if c in s:
            s[c] = pd.to_numeric(s[c], errors="coerce")
    med = float(s["reads_after_trim"].median()) if "reads_after_trim" in s else float(s["reads_in"].median())
    rows = []
    for _, r in s.iterrows():
        fp = q.fastp.get(str(r["sample_id"]), {})
        after = fp.get("after_filtering", {})
        kept = r.get("reads_after_trim")
        row = {"sample_id": r["sample_id"], "layout": r.get("layout", ""), "reads_in": r["reads_in"],
               "reads_after_trim": kept, "retained": kept / r["reads_in"] if kept and r["reads_in"] else None,
               "mapping_rate": r["mapping_rate"], "q30": after.get("q30_rate"), "gc": after.get("gc_content"),
               "flags": []}
        if row["mapping_rate"] is not None and row["mapping_rate"] < LOW_MAPPING:
            row["flags"].append(f"mapping rate < {LOW_MAPPING:.0%}")
        depth = kept if kept else r["reads_in"]
        if depth < LOW_READS:
            row["flags"].append(f"< {LOW_READS / 1e6:g} M reads")
        if depth * OUTLIER_FOLD < med:
            row["flags"].append(f"{OUTLIER_FOLD:g}x below the median depth")
        rows.append(row)
    mr = s["mapping_rate"]
    obs = [f"{len(rows)} libraries; {num(s['reads_in'].median() / 1e6, 1)} M reads per library (median), "
           f"mapping rate {pct(mr.min())} to {pct(mr.max())} (median {pct(mr.median())})."]
    if "reads_after_trim" in s:
        ret = (s["reads_after_trim"] / s["reads_in"]).median()
        obs.append(f"Trimming kept {pct(ret)} of reads (median).")
    flagged = [r for r in rows if r["flags"]]
    if flagged:
        obs.append(f"{len(flagged)} librar{'y' if len(flagged) == 1 else 'ies'} flagged: "
                   + "; ".join(f"{r['sample_id']} ({', '.join(r['flags'])})" for r in flagged) + ".")
    caveats = [f"Flags are rules of thumb (mapping rate < {LOW_MAPPING:.0%}, < {LOW_READS / 1e6:g} M reads); a "
               "flagged library may still be fine for gene-level counts. Check whether it stands apart in the PCA "
               "before removing it."]
    return {"source": "quant", "rows": rows, "observations": obs, "caveats": caveats, "n_flagged": len(flagged),
            "has_q30": any(r["q30"] is not None for r in rows)}


# ---- differential expression -------------------------------------------------------------------

def _r2(values: pd.Series, groups: pd.Series) -> float:
    """Share of the variance of values explained by the groups (one-way ANOVA R^2)."""
    total = float(((values - values.mean()) ** 2).sum())
    if total == 0:
        return 0.0
    within = float(sum(((g - g.mean()) ** 2).sum() for _, g in values.groupby(groups)))
    return 1 - within / total


def de_run(run: DeRun, fdr: float, max_genes: int) -> dict[str, Any]:
    s, d = run.summary, run.design
    terms = d.get("terms") or []
    out: dict[str, Any] = {"label": run.label, "formula": d.get("formula"), "terms": terms,
                           "levels": d.get("levels", {}), "reference_levels": d.get("reference_levels", {}),
                           "samples": d.get("samples", []), "left_out": d.get("samples_left_out", []),
                           "genes": d.get("genes", {}), "organism": s.get("organism"), "layout": s.get("layout"),
                           "notes": s.get("notes", []), "pca_png": run.files.png(run.file("pca.png")), "contrasts": []}
    obs, caveats = [], []
    n = len(out["samples"])
    per_level = {t: pd.Series([x[t] for x in out["samples"]]).value_counts().to_dict() for t in terms}
    out["per_level"] = per_level
    smallest = min((min(v.values()) for v in per_level.values() if v), default=0)
    obs.append(f"{n} samples, design {out['formula']}: " + "; ".join(
        f"{t} " + ", ".join(f"{k} ({v})" for k, v in per_level[t].items()) for t in terms) + ".")
    g = out["genes"]
    if g:
        obs.append(f"{num(g.get('tested'))} of {num(g.get('input'))} genes tested (>= {g.get('min_count')} counts in "
                   f">= {g.get('min_samples')} samples).")
    if run.pca is not None and "PC2" in run.pca:
        var = (s.get("pca") or {}).get("variance_explained") or []
        for i, pc in enumerate(("PC1", "PC2")):
            if pc not in run.pca or i >= len(var):
                continue
            r2 = {t: _r2(run.pca[pc], run.pca[t].astype(str)) for t in terms if t in run.pca}
            if r2:
                best = max(r2, key=r2.get)
                obs.append(f"{pc} ({pct(var[i], 0)} of the variance of the 500 most variable genes) "
                           + (f"separates the samples by {best} (R² {r2[best]:.2f})." if r2[best] >= 0.5 else
                              f"is not explained by a design term (best: {best}, R² {r2[best]:.2f})."))
        out["pca_variance"] = var
    if smallest and smallest < 3:
        caveats.append(f"A group has {smallest} sample(s): DESeq2 can test it, but the dispersion rests on very "
                       "few replicates and only large, consistent changes are reliable.")
    elif smallest:
        caveats.append(f"The smallest group has {smallest} samples: modest effects can be missed; a gene that is "
                       "not significant is not shown to be unchanged.")
    caveats.append("Fold changes are between groups of libraries. Pooled libraries (several animals per library) "
                   "average out differences between the animals, and samples processed in batches can differ "
                   "for technical reasons: include the batch in the design when it is known.")

    for c in s.get("contrasts", []):
        res = run.files.csv(run.file(f"{c['name']}.csv"), "gene_id", "padj", keep_default_na=False,
                       na_values=[""], dtype={"gene_id": str, "gene_name": str})
        e = {**c, "volcano": run.files.png(run.file(f"{c['name']}_volcano.png")),
             "ma": run.files.png(run.file(f"{c['name']}_ma.png")),
             "up": [], "down": [], "gsea": []}
        if res is not None:
            res["gene_name"] = res["gene_name"].fillna("")
            lfc = "log2_fold_change_shrunk" if "log2_fold_change_shrunk" in res else "log2_fold_change"
            sig = res[res["padj"] < fdr]
            e["n_significant"] = len(sig)
            e["n_up"] = int((sig["log2_fold_change"] > 0).sum())
            e["n_down"] = int((sig["log2_fold_change"] < 0).sum())
            e["lfc_column"] = lfc

            def rows(df):
                return [{"gene_id": r["gene_id"], "gene_name": r["gene_name"], "base_mean": r["base_mean"],
                         "lfc": r[lfc], "padj": r["padj"]} for _, r in df.iterrows()]
            e["up"] = rows(sig[sig["log2_fold_change"] > 0].nsmallest(max_genes, "padj"))
            e["down"] = rows(sig[sig["log2_fold_change"] < 0].nsmallest(max_genes, "padj"))
        gs = run.files.csv(run.file(f"{c['name']}_gsea.csv"), "collection", "term", "nes", "fdr")
        if gs is not None:
            for coll, sub in gs.groupby("collection", sort=False):
                sig_sets = sub[sub["fdr"] < fdr]
                e["gsea"].append({
                    "collection": coll, "n_sets": len(sub), "n_significant": len(sig_sets),
                    "n_up": int((sig_sets["nes"] > 0).sum()), "n_down": int((sig_sets["nes"] < 0).sum()),
                    "up": sub[sub["nes"] > 0].nsmallest(max_genes // 2 or 1, "fdr").to_dict("records"),
                    "down": sub[sub["nes"] < 0].nsmallest(max_genes // 2 or 1, "fdr").to_dict("records")})
        txt = (f"{c['column']} {c['level']} vs {c['reference']} ({c['n_samples_level']} vs {c['n_samples_reference']} "
               f"samples): {e.get('n_up', 0)} genes up and {e.get('n_down', 0)} down at FDR {fdr:g}")
        if e["up"] or e["down"]:
            tops = [x["gene_name"] or x["gene_id"] for x in (e["up"][:3] + e["down"][:3])]
            txt += "; strongest: " + ", ".join(tops)
        e["observation"] = txt + "."
        for gse in e["gsea"]:
            if gse["n_significant"]:
                lead = [x["term"] for x in (gse["up"][:2] + gse["down"][:2]) if x["fdr"] < fdr]
                e.setdefault("gsea_observations", []).append(
                    f"{gse['collection']}: {gse['n_up']} sets enriched among up-regulated genes, {gse['n_down']} "
                    f"among down-regulated (FDR {fdr:g})" + (f", e.g. {', '.join(lead)}" if lead else "") + ".")
        out["contrasts"].append(e)
        obs.append(e["observation"])
    if any(not c.get("lfc_shrunk") for c in s.get("contrasts", [])):
        caveats.append("Some contrasts have no shrunk fold changes (the reference is not the factor's base level); "
                       "their unshrunk estimates are noisy for genes with few counts.")
    if any(e["gsea"] for e in out["contrasts"]):
        caveats.append("GSEA tests whether a set's genes sit at the top or bottom of the ranking by Wald statistic; "
                       "an enriched set is evidence about many modest changes together, not about any single gene, "
                       "and overlapping GO sets report the same genes several times.")
    out["observations"], out["caveats"] = obs, caveats
    return out


# ---- methods, citations, credits ---------------------------------------------------------------

# Printed first in every report's citation section, above the methods (the paper is "pennsieve" in references.json).
PLATFORM_STATEMENT = (
    "Citing Pennsieve. If you used Pennsieve to store, analyse, share or publish data reported in a publication, "
    "please cite: Goldblum Z, et al. Pennsieve: A Collaborative Platform for Translational Neuroscience and Beyond. "
    "Scientific Data 2025;12:1834. doi:10.1038/s41597-025-06075-5. Please also cite the datasets you used (their "
    "DOIs) and the methods listed below.")


def platform_citation() -> dict[str, Any]:
    r = _references()["pennsieve"]
    return {"statement": PLATFORM_STATEMENT, "doi": r["doi"], "pmid": r.get("pmid") or "", "text": _ref_text(r)}


def _references() -> dict[str, dict[str, Any]]:
    return json.loads(REFERENCES.read_text(encoding="utf-8"))


def _norm_doi(doi: str) -> str:
    return re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", doi.strip(), flags=re.I).lower()


def _ref_text(r: dict[str, Any]) -> str:
    where = str(r.get("volume") or "") + (f"({r['issue']})" if r.get("issue") else "") + \
        (f":{r['pages']}" if r.get("pages") else "")
    return f"{r['short_authors']}. {r['title']}. {r['journal_short']} {r['year']}" + (f";{where}" if where else "") + "."


def methods(inp: Inputs, own: dict[str, Any]) -> dict[str, Any]:
    """Steps from every provenance.json, and the citations they list (objects with doi/text, or DOI strings)."""
    known = {_norm_doi(r["doi"]): r for r in _references().values()}
    refs: list[dict[str, Any]] = []
    index: dict[str, int] = {}

    def add(c) -> int | None:
        if isinstance(c, str):
            c = {"doi": c}
        if not isinstance(c, dict) or not (c.get("doi") or c.get("text")):
            return None
        doi = _norm_doi(str(c.get("doi") or ""))
        k = doi or str(c.get("text"))
        if k not in index:
            r = known.get(doi)
            # The verified record wins over a step's own text; anything else is shown as given, marked.
            text = _ref_text(r) if r else str(c.get("text") or f"doi:{doi}")
            refs.append({"doi": doi, "pmid": (r or {}).get("pmid") or str(c.get("pmid") or ""), "text": text,
                         "verified": r is not None})
            index[k] = len(refs)
        return index[k]

    steps = []
    provs = list(inp.provenances) + [(None, own)]
    for path, p in provs:
        nums = sorted({n for n in (add(c) for c in p.get("citations") or []) if n})
        params = p.get("parameters") or {}
        steps.append({"app": p.get("app"), "version": p.get("version"),
                      "where": inp.rel(path.parent) if path is not None else "this report",
                      "tools": ", ".join(f"{k} {v}" for k, v in (p.get("tools") or {}).items() if v),
                      "parameters": ", ".join(f"{k}={v}" for k, v in params.items() if v not in ("", None)),
                      "citations": nums})
    return {"steps": steps, "references": refs}


DATASETS = [
    # Discover 427: its pool-<id> folder names, as bulk-de records them (generic names like N1, T1 do not count)
    {"match": {"pool-FN1", "pool-FN2", "pool-FN3", "pool-FT1", "pool-FT2", "pool-FT3", "pool-N1", "pool-N2",
               "pool-N3", "pool-T1", "pool-T2", "pool-MT"}, "ref": "tg_data",
     "paper": "tg_paper", "licence": "CC BY 4.0", "url": "https://sparc.science/datasets/427"},
]


def credits(inp: Inputs, data_credit: str) -> list[dict[str, Any]]:
    out = []
    if data_credit:
        out.append({"text": data_credit, "doi": "", "licence": "", "paper": None})
    samples = {str(x) for r in inp.de_runs for x in r.de.get("samples") or []}
    refs = _references()
    for d in DATASETS:
        if len(samples) >= 6 and samples <= d["match"]:
            r = refs[d["ref"]]
            authors = ", ".join(f"{a['family']} {''.join(x[0] for x in a['given'].split())}" for a in r["authors"])
            out.append({"text": f"{authors}. {r['title']}. {r['journal']} ({r['year']}).", "doi": r["doi"],
                        "licence": d["licence"], "url": d["url"],
                        "paper": _ref_text(refs[d["paper"]]) + f" doi:{refs[d['paper']]['doi']}"})
    return out
