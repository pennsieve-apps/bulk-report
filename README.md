# bulk-report

The last step of the bulk RNA-seq workflows: one HTML page that a biologist can open offline
and a reviewer can check, from bulk-de's differential expression and, when linked,
rnaseq-quant's sample QC.

```
A. rnaseq-quant (+ index) -> bulk-de -> bulk-report      (link rnaseq-quant's counts to `counts` too)
B.        deposited counts -> bulk-de -> bulk-report
```

The page has, in order: the **design** (formula, samples per level, reference levels, samples
left out and why); **sample QC** (reads, trimming, salmon mapping rate, fastp Q30/GC per library,
with flags; or, without rnaseq-quant, counts and genes detected per library); per comparison the
**PCA** with a sentence on which design term PC1/PC2 follow, and per contrast the **volcano and
MA plots**, the **top genes** up and down, and the **enriched gene sets** per collection; then
**methods** (every step's app, version, parameters and tools, from its `provenance.json`),
**citations** gathered from every step's provenance and numbered, and **data credits**. Each
section states observations backed by numbers and the caveats that go with them (few replicates,
pooled libraries, additive designs, GSEA on overlapping sets).

The page is a single file: figures are embedded as PNG data URIs, there are no scripts, fonts or
external links other than DOIs. It is generated with Jinja2's autoescaping, so text from the
inputs (sample names, sheet values) cannot inject markup.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `de` | input (`placement: subfolder`, required) | bulk-de `<LABEL>/` folders | Several links are fine: every folder with `design.json` + `summary.json` is one comparison. |
| `counts` | input (`placement: subfolder`, optional) | rnaseq-quant's `samples.csv` (`reads_in`, `reads_after_trim`, `mapping_rate`), `qc/<sample>.fastp.json`, `provenance.json` | Read-level QC. |
| `report` | output | `index.html`, `summary.json`, `citations.csv`, `provenance.json` | |

Without links, everything is found by name anywhere under `INPUT_DIR` (symlinks followed; a
link below an entry may not leave that entry's tree). Every `provenance.json` with an `app`
under the input contributes a methods row and its `citations`.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `REPORT_TITLE` | string | `Bulk RNA-seq Report` | Page title. |
| `FDR` | number | `0.05` | What the page counts and lists as significant, applied to bulk-de's adjusted p-values (display only). |
| `MAX_GENES` | integer | `15` | Rows per up/down gene table; gene sets per direction: half. |
| `DATA_CREDIT` | string | *(empty)* | Credit line for the input data (authors, title, repository, DOI, licence). Discover 427 is recognised from its library names without it. |

Exit codes: `0` success; `1` no bulk-de output and no rnaseq-quant `samples.csv`; `2` bad
parameters.

## Outputs

- `index.html`: the page.
- `summary.json`: title, FDR; sample QC (source, observations, caveats, flagged count); per
  comparison the formula, samples, observations, caveats and per contrast n up/down, top genes,
  gene-set counts; the steps with their citation numbers; the numbered citations
  (`doi`, `pmid`, `text`); data credits; section errors.
- `citations.csv`: `number, citation, doi, pmid`.
- `provenance.json`: app, version, tools, parameters, inputs, started, finished, samples (QC
  flags), peak RSS, citations.

**Citations.** Each step's `provenance.json` lists `citations` as `{doi, pmid, text}` (the
single-cell and bulk conventions); bare DOI strings are accepted too. They are deduplicated by
DOI (case-insensitive, `doi:`/`https://doi.org/` prefixes removed) and numbered in pipeline
order. An entry without text is completed from `processor/references.json` (bulk-de's and
rnaseq-quant's methods, verified; see Citations) or shown as its DOI.

**Sample QC flags** (rules of thumb, stated on the page): mapping rate < 50 %, < 5 M reads after
trimming, depth 3x below the median; without rnaseq-quant, a library 3x below the median count.

## Running it locally

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python     # pytest on synthetic upstream outputs, linked and flat layouts
make run PYTHON=.venv/bin/python      # image (amd64) over the fixture -> data/output/index.html
```

The tests check the page is self-contained (no scripts, no external `src`/`href`), every figure
embedded, the design / QC / contrast / GSEA content, the QC flag, citations gathered from two
provenance files (objects and bare DOIs, one unknown DOI) and deduplicated, the Discover 427
credit, `DATA_CREDIT`, escaping of hostile text, bad parameters, no input, and the Lambda handler.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 2048 MB, 5 min timeout. It reads summary files and,
per contrast, the gene table for the top genes; nothing per read. Measured with
`docker run --platform linux/amd64 --memory 2g --cpus 1 --network none --read-only --tmpfs /tmp`
(amd64 emulation on an M-series Mac): Discover 427 (bulk-de with 2 contrasts x 13 k genes and
GSEA on 3,769 sets, plus rnaseq-quant QC for 12 libraries): **2.0 s, 156 MB peak RSS**, page
367 KB. Image: 328 MB unpacked, 109 MB compressed.

## Results on the demo data

`data/work/tg-tumour/report/report/index.html` in the workflow repository, from
`data/work/tg-tumour/de/de/` (bulk-de on the authors' counts, `~ sex + condition`) and
`data/work/tg-tumour/quant/counts/` (rnaseq-quant on 1 M-read subsamples). It reports 12
libraries with 87 % median mapping rate (64-89 %), all flagged "< 5 M reads" because of the
subsampling; PC1 (68 %) not explained by sex or condition; tumour vs naïve: no gene at FDR 0.05,
2 Hallmark sets (TNF-α/NF-κB, allograft rejection) and 23 GO:BP sets (MHC class II antigen
presentation) up; male vs female: 30 up / 652 down, with Y genes, Xist and the female-higher
immune markers (C1qa/b/c, Lyz2, Fcer1g, Tyrobp) the paper describes; 16 citations from the
three steps; the Discover 427 credit with its describing paper.

## Limitations

- The page describes what the upstream steps computed; it does not re-test anything. `FDR`
  changes what is counted, not the tests.
- Observations are templated from numbers (which term PC1 follows, how many genes and sets); they
  are not an interpretation of the biology.
- Data credits are automatic only for Discover 427; otherwise set `DATA_CREDIT`.
- Citations without text and with a DOI it does not know are shown as the bare DOI.

## Licences

Our code: Apache-2.0 (`LICENSE`). The image holds only permissive components (Jinja2 BSD-3,
pandas/NumPy BSD-3, Matplotlib PSF-based); see `THIRD_PARTY.md`.

## Credits

- **Matplotlib** (John D. Hunter and the Matplotlib developers): the QC figure.
- **Jinja2** (Pallets): the HTML template.
- The methods the page reports are credited by the steps that ran them (bulk-de: PyDESeq2,
  DESeq2, apeglm, GSEA, GSEApy, MSigDB, GO; rnaseq-quant: salmon, fastp, tximport, GENCODE),
  and the page lists them.
- **Demo data:** Discover dataset 427 by Alexei V. Tumanov, Sergey A. Shein, Zhao Lai, Yidong
  Chen and Shivani B. Ruparel (RE-JOIN; CC BY 4.0), described by Ibrahim et al. 2023.

## Citations

Checked on 2026-10-31 against Crossref and PubMed (esearch on the DOI, then esummary), the
dataset against DataCite and the Discover API; `scripts/verify_citations.py` rewrites
`processor/references.json`, which also holds the records used to complete DOI-only citations
(salmon, fastp, tximport, GENCODE 2023 and bulk-de's methods).

- Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput Sci Eng* 9(3):90-95 (2007).
  doi:10.1109/MCSE.2007.55. (Crossref; not in PubMed)

Demo data: Tumanov AV, Shein SA, Lai Z, Chen Y, Ruparel SB. Sex-dependent differences in the
genomic profile of lingual sensory neurons in naive and tongue-tumor bearing mice. SPARC Portal
(2025). doi:10.26275/gaap-zlpj. CC BY 4.0. Described by Ibrahim T, Wu P, Wang LJ, et al.
*Sci Rep* 13(1):13117 (2023). doi:10.1038/s41598-023-40380-6. PMID 37573456.
