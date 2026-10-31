"""
Synthetic upstream outputs in the shape bulk-de and rnaseq-quant write them:

  linked/de/condition/   a bulk-de run (~ sex + condition, contrasts condition:T:N and sex:M:F)
  linked/de/sex-only/    a second bulk-de run (~ sex), to check several runs in one port
  linked/counts/         rnaseq-quant's samples.csv (reads, trimming, mapping rate), qc/*.fastp.json and a
                         provenance.json whose citations are bare DOI strings
  flat/                  the same, everything side by side without port folders

8 libraries named as Discover 427's (FN1, FT1, ..., so the data credit is recognised); FT2 has a
mapping rate of 31 % (flagged).

    python -m tests.make_fixture data/fixture
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = [("FN1", "F", "N"), ("FN2", "F", "N"), ("FT1", "F", "T"), ("FT2", "F", "T"),
           ("N1", "M", "N"), ("N2", "M", "N"), ("T1", "M", "T"), ("T2", "M", "T")]
DE_CITES = [{"doi": "10.1093/bioinformatics/btad547", "pmid": "37669147",
             "text": "Muzellec B, Teleńczuk M, Cabeli V, Andreux M. PyDESeq2: a python package for bulk RNA-seq "
                     "differential expression analysis. Bioinformatics 2023;39(9):btad547."},
            {"doi": "10.1186/s13059-014-0550-8", "pmid": "25516281",
             "text": "Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data "
                     "with DESeq2. Genome Biol 2014;15(12):550."},
            {"doi": "10.1073/pnas.0506580102", "pmid": "16199517",
             "text": "Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis. PNAS 2005."}]


def _png(path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(2, 1.5))
    ax.plot([0, 1], [0, 1])
    fig.savefig(path, dpi=50)
    plt.close(fig)


def de_run(out: Path, label: str, terms: list[str], contrasts: list[tuple[str, str, str]], rng) -> None:
    out.mkdir(parents=True)
    samples = [{"sample_id": f"pool-{s}", "sex": sex, "condition": cond} for s, sex, cond in SAMPLES]
    design = {"formula": "~ " + " + ".join(terms), "terms": terms,
              "contrasts": [":".join(c) for c in contrasts],
              "reference_levels": {"sex": "F", "condition": "N"}, "levels": {"sex": ["F", "M"], "condition": ["N", "T"]},
              "samples": [{k: s[k] for k in ["sample_id"] + terms} for s in samples], "samples_left_out": [],
              "genes": {"input": 300, "tested": 250, "min_count": 10, "min_samples": 4}}
    (out / "design.json").write_text(json.dumps(design))
    pd.DataFrame([{**s, "library_size": int(rng.integers(8e6, 2e7)), "n_genes_detected": 14000,
                   "size_factor": 1.0, "used": True} for s in samples]).to_csv(out / "samples.csv", index=False)
    pc1 = [3.0 if s["condition"] == "T" else -3.0 for s in samples]
    pd.DataFrame([{**s, "PC1": p + rng.normal(0, 0.3), "PC2": rng.normal()} for s, p in zip(samples, pc1)]).to_csv(
        out / "pca.csv", index=False)
    _png(out / "pca.png")
    entries = []
    for col, lv, ref in contrasts:
        name = f"{col}-{lv}-vs-{ref}"
        genes = pd.DataFrame({"gene_id": [f"ENSMUSG{i:011d}" for i in range(250)],
                              "gene_name": [f"Gene{i}" for i in range(250)],
                              "base_mean": rng.uniform(10, 1000, 250), "log2_fold_change": rng.normal(0, 0.3, 250),
                              "lfc_se": 0.2, "stat": rng.normal(0, 1, 250), "pvalue": rng.uniform(0.01, 1, 250),
                              "padj": rng.uniform(0.2, 1, 250)})
        genes.loc[:9, "log2_fold_change"] = 2.0
        genes.loc[10:14, "log2_fold_change"] = -1.5
        genes.loc[:14, "padj"] = np.linspace(1e-12, 1e-3, 15)
        genes["log2_fold_change_shrunk"] = genes["log2_fold_change"] * 0.8
        genes.to_csv(out / f"{name}.csv", index=False)
        _png(out / f"{name}_volcano.png")
        _png(out / f"{name}_ma.png")
        pd.DataFrame({"collection": ["hallmark"] * 3, "term": ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_MYC",
                                                                "HALLMARK_APOPTOSIS"],
                      "size": [150, 180, 160], "es": [0.6, -0.5, 0.1], "nes": [2.1, -1.8, 0.4],
                      "pvalue": [0.001, 0.002, 0.8], "fdr": [0.004, 0.01, 0.9], "n_leading_edge": [30, 25, 10],
                      "leading_edge": ["Il1b;Ccl2;Cxcl1", "Myc;Npm1", "Bax"]}).to_csv(out / f"{name}_gsea.csv",
                                                                                       index=False)
        entries.append({"contrast": f"{col}:{lv}:{ref}", "name": name, "column": col, "level": lv, "reference": ref,
                        "n_samples_level": 4, "n_samples_reference": 4, "n_up": 10, "n_down": 5,
                        "lfc_shrunk": True})
    (out / "summary.json").write_text(json.dumps({"label": label, "formula": design["formula"], "fdr": 0.05,
                                                  "organism": "mouse", "layout": "two-column",
                                                  "pca": {"variance_explained": [0.42, 0.18]},
                                                  "contrasts": entries, "notes": []}))
    (out / "de.json").write_text(json.dumps({"layout": "two-column", "samples": [s["sample_id"] for s in samples]}))
    (out / "provenance.json").write_text(json.dumps({"app": "bulk-de", "version": "0.1.0", "tools": {"pydeseq2": "0.5.4"},
                                                     "parameters": {"DESIGN": design["formula"]},
                                                     "citations": DE_CITES}))


def quant(out: Path, rng) -> None:
    (out / "qc").mkdir(parents=True)
    rows = []
    for s, sex, cond in SAMPLES:
        reads = int(rng.integers(15e6, 25e6))
        rows.append({"sample_id": s, "files": f"pool-{s}/{s}_R1.fastq.gz", "layout": "single", "reads_in": reads,
                     "reads_after_trim": int(reads * 0.97), "mapping_rate": 0.31 if s == "FT2" else 0.85,
                     "library_type": "SR", "sex": sex, "condition": cond})
        (out / "qc" / f"{s}.fastp.json").write_text(json.dumps(
            {"summary": {"after_filtering": {"q30_rate": 0.93, "gc_content": 0.49}}}))
    pd.DataFrame(rows).to_csv(out / "samples.csv", index=False)
    (out / "provenance.json").write_text(json.dumps({"app": "rnaseq-quant", "version": "0.1.0",
                                                     "tools": {"salmon": "1.10.3"}, "parameters": {},
                                                     "citations": ["10.1038/nmeth.4197", "doi:10.1093/bioinformatics/bty560", "10.1002/imt2.70078",
                                                                   "10.1186/s13059-014-0550-8", "10.9999/unknown"]}))


def make_fixture(root: Path) -> Path:
    rng = np.random.default_rng(0)
    linked = root / "linked"
    de_run(linked / "de" / "condition", "condition", ["sex", "condition"],
           [("condition", "T", "N"), ("sex", "M", "F")], rng)
    de_run(linked / "de" / "sex-only", "sex-only", ["sex"], [("sex", "M", "F")], rng)
    quant(linked / "counts", rng)
    flat = root / "flat"
    for src, name in ((linked / "de" / "condition", "a"), (linked / "de" / "sex-only", "b"), (linked / "counts", "c")):
        shutil.copytree(src, flat / name)
    return root


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    make_fixture(Path(ap.parse_args().out))
