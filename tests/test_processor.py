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
    assert not re.search(r'(src|href)="(?!data:|#|https://doi\.org/)', html)
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
    assert salmon["text"].startswith("Patro R") and salmon["pmid"] == "28263959"
    assert next(x for x in cites if x["doi"] == "10.9999/unknown")["text"] == "doi:10.9999/unknown"
    assert {s["app"] for s in summary["steps"]} == {"bulk-de", "rnaseq-quant", "bulk-report"}
    # data credit for Discover 427 recognised from the library names
    assert summary["data"] and summary["data"][0]["doi"] == "10.26275/gaap-zlpj"
    assert "Tumanov AV" in summary["data"][0]["text"] and "Ruparel SB" in summary["data"][0]["text"]
    assert (tmp_path / "citations.csv").read_text().count("\n") == len(cites) + 1
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
