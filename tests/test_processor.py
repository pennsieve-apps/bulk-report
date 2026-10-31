"""End-to-end tests on synthetic upstream outputs (tests/make_fixture.py), no Docker."""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from tests.make_fixture import make_fixture

PARAMS = ("REPORT_TITLE", "FDR", "MAX_GENES", "DATA_CREDIT")


@pytest.fixture(scope="session")
def fx(tmp_path_factory):
    return make_fixture(tmp_path_factory.mktemp("fixture"))


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for k, v in params.items():
        monkeypatch.setenv(k, str(v))
    from processor.main import main

    return main()


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_report(fx, tmp_path, monkeypatch, layout):
    assert run(fx / layout, tmp_path, monkeypatch) == 0
    html = (tmp_path / "index.html").read_text()
    # self-contained: no scripts, no external resources; figures embedded
    assert "<script" not in html.lower()
    assert not re.search(r'(src|href)="(?!data:|#|https://doi\.org/|https://sparc\.science/datasets/)', html)
    assert "Content-Security-Policy" in html and "default-src 'none'" in html
    assert html.count('src="data:image/png;base64,') >= 1 + 1 + 3 * 2   # QC, 2 PCA, 3 contrasts x (volcano, MA)
    for text in ("~ sex + condition", "condition: T vs N", "sex: M vs F", "HALLMARK_INFLAMMATORY_RESPONSE",
                 "Gene0", "mapping rate &lt; 50%", "FT2"):
        assert text in html, text
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert [r["label"] for r in summary["differential_expression"]] == ["condition", "sex-only"]
    c = summary["differential_expression"][0]["contrasts"][0]
    assert c["n_up"] == 10 and c["n_down"] == 5 and c["top_up"][0] == "Gene0"
    assert summary["sample_qc"]["source"] == "quant" and summary["sample_qc"]["n_flagged"] == 1
    assert any("PC1" in o and "condition" in o for o in summary["differential_expression"][0]["observations"])
    # citations: gathered from every provenance.json, deduplicated, DOI strings resolved from references.json
    cites = summary["citations"]
    dois = [x["doi"] for x in cites]
    assert len(dois) == len(set(dois))
    assert {"10.1093/bioinformatics/btad547", "10.1038/nmeth.4197", "10.1093/bioinformatics/bty560",
            "10.9999/unknown", "10.1109/mcse.2007.55"} <= set(dois)
    salmon = next(x for x in cites if x["doi"] == "10.1038/nmeth.4197")
    assert salmon["text"].startswith("Patro R") and salmon["pmid"] == "28263959" and salmon["verified"]
    unknown = next(x for x in cites if x["doi"] == "10.9999/unknown")
    assert unknown["text"] == "doi:10.9999/unknown" and not unknown["verified"]
    # the verified record wins over a step's own (shorter) text
    gsea = next(x for x in cites if x["doi"] == "10.1073/pnas.0506580102")
    assert gsea["text"].startswith("Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: A")
    assert "not checked by this report" in html
    statement = ("Citing Pennsieve. If you used Pennsieve to store, analyse, share or publish data reported in a "
                 "publication, please cite: Goldblum Z, et al. Pennsieve: A Collaborative Platform for Translational "
                 "Neuroscience and Beyond. Scientific Data 2025;12:1834. doi:10.1038/s41597-025-06075-5. Please also "
                 "cite the datasets you used (their DOIs) and the methods listed below.")
    cit = html[html.index('id="citations"'):]
    assert statement in cit and cit.index(statement) < cit.index("<ol")
    assert summary["platform_citation"]["statement"] == statement
    assert summary["platform_citation"]["doi"] == "10.1038/s41597-025-06075-5"
    assert {s["app"] for s in summary["steps"]} == {"bulk-de", "rnaseq-quant", "bulk-report"}
    # data credit for Discover 427 recognised from the library names
    assert summary["data"] and summary["data"][0]["doi"] == "10.26275/gaap-zlpj"
    assert "Tumanov AV" in summary["data"][0]["text"] and "Ruparel SB" in summary["data"][0]["text"]
    import pandas as pd

    table = pd.read_csv(tmp_path / "citations.csv", dtype=str, keep_default_na=False)
    assert len(table) == len(cites) + 1 and table.at[0, "kind"] == "platform"
    assert table.at[0, "citation"].startswith("Goldblum Z, Xu Z, Shi H, et al. Pennsieve")
    assert (table["kind"].iloc[1:] == "method").all()
    fastp1 = next(x for x in cites if x["doi"] == "10.1002/imt2.70078")
    assert fastp1["verified"] and fastp1["text"].startswith("Chen S. fastp 1.0") and fastp1["pmid"] == "41112039"
    prov = json.loads((tmp_path / "provenance.json").read_text())
    for k in ("app", "version", "tools", "parameters", "inputs", "started", "finished", "samples", "peak_rss_mb"):
        assert k in prov
    assert prov["app"] == "bulk-report" and prov["citations"][0]["text"].startswith("Hunter JD")


def test_de_only_and_credit(fx, tmp_path, monkeypatch):
    """bulk-de alone (no rnaseq-quant): QC from library sizes; DATA_CREDIT shown."""
    assert run(fx / "linked" / "de", tmp_path, monkeypatch, DATA_CREDIT="Someone A. Some data. doi:10.1/x",
               MAX_GENES=3) == 0
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["sample_qc"]["source"] == "counts"
    assert summary["data"][0]["text"].startswith("Someone A")
    assert len(summary["differential_expression"][0]["contrasts"][0]["top_up"]) == 3


def test_hostile_text_is_escaped(fx, tmp_path, monkeypatch):
    assert run(fx / "linked", tmp_path, monkeypatch, REPORT_TITLE="<script>alert(1)</script>") == 0
    assert "<script>" not in (tmp_path / "index.html").read_text()


def test_bad_parameters_and_missing_inputs(fx, tmp_path, monkeypatch):
    assert run(fx / "linked", tmp_path, monkeypatch, FDR="0") == 2
    assert run(fx / "linked", tmp_path, monkeypatch, MAX_GENES="x") == 2
    (tmp_path / "empty").mkdir()
    assert run(tmp_path / "empty", tmp_path / "out", monkeypatch) == 1


def test_handler(fx, tmp_path, monkeypatch):
    for k in PARAMS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("REPORT_TITLE", "stale")
    from processor.handler import handler

    assert handler({"inputDir": str(fx / "linked"), "outputDir": str(tmp_path)}, None)["status"] == "ok"
    assert json.loads((tmp_path / "summary.json").read_text())["title"] == "Bulk RNA-seq Report"
    with pytest.raises(RuntimeError):
        handler({"inputDir": str(fx / "linked"), "outputDir": str(tmp_path), "FDR": 5}, None)


def _copy_run(src: Path, dst: Path) -> Path:
    import shutil

    shutil.copytree(src, dst)
    return dst


def test_platform_symlinks(fx, tmp_path, monkeypatch):
    """INPUT_DIR/de/ is a real folder whose entries are symlinks to the upstream outputs elsewhere."""
    store = tmp_path / "store"
    _copy_run(fx / "linked" / "de" / "condition", store / "run1" / "condition")
    _copy_run(fx / "linked" / "counts", store / "quant")
    inp = tmp_path / "in"
    (inp / "de").mkdir(parents=True)
    (inp / "de" / "condition").symlink_to(store / "run1" / "condition")
    (inp / "counts").symlink_to(store / "quant")
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert [r["label"] for r in summary["differential_expression"]] == ["condition"]
    assert summary["sample_qc"]["source"] == "quant"


def test_files_next_to_design_json_cannot_escape(fx, tmp_path, monkeypatch):
    """A link inside an upstream folder that points out of the input is neither read nor embedded."""
    import shutil

    inp = tmp_path / "in"
    run_dir = _copy_run(fx / "linked" / "de" / "condition", inp / "de" / "condition")
    secret = tmp_path / "secret"
    secret.mkdir()
    shutil.copy(run_dir / "pca.png", secret / "secret.png")
    (secret / "secret.csv").write_text("sample_id,library_size,n_genes_detected,size_factor,used\nSECRET,1,1,1,True\n")
    (run_dir / "pca.png").unlink()
    (run_dir / "pca.png").symlink_to(secret / "secret.png")
    (run_dir / "samples.csv").unlink()
    (run_dir / "samples.csv").symlink_to(secret / "secret.csv")
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    html = (tmp_path / "out" / "index.html").read_text()
    assert "SECRET" not in html
    assert html.count('src="data:image/png') == 4          # 2 contrasts x (volcano, MA); no PCA, no QC figure


def test_one_bad_section_does_not_sink_the_page(fx, tmp_path, monkeypatch):
    inp = tmp_path / "in"
    run_dir = _copy_run(fx / "linked" / "de" / "condition", inp / "de" / "condition")
    _copy_run(fx / "linked" / "de" / "sex-only", inp / "de" / "sex-only")
    s = json.loads((run_dir / "summary.json").read_text())
    s["contrasts"] = ["not a contrast"]
    (run_dir / "summary.json").write_text(json.dumps(s))
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    summary = json.loads((tmp_path / "out" / "summary.json").read_text())
    assert [r["label"] for r in summary["differential_expression"]] == ["sex-only"]
    assert summary["errors"] and "condition" in summary["errors"][0]["section"]
    assert "could not be rendered" in (tmp_path / "out" / "index.html").read_text()


def test_generic_library_names_get_no_credit(fx, tmp_path, monkeypatch):
    """Library names like N1, T1 alone must not credit Discover 427."""
    inp = tmp_path / "in"
    run_dir = _copy_run(fx / "linked" / "de" / "condition", inp / "de" / "condition")
    de = json.loads((run_dir / "de.json").read_text())
    de["samples"] = [x.removeprefix("pool-") for x in de["samples"]]
    (run_dir / "de.json").write_text(json.dumps(de))
    assert run(inp, tmp_path / "out", monkeypatch) == 0
    assert json.loads((tmp_path / "out" / "summary.json").read_text())["data"] == []


def test_gsea_p_zero_shown_as_bound():
    from processor.sections import gsea_p

    assert gsea_p(0.0) == "< 1/1001" and gsea_p(0.004) == "0.004"
