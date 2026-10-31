# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image contains the following third-party components,
installed from PyPI at the pinned versions in `requirements.txt`; their licences apply to the
image. Nothing here is GPL. No data files are bundled except `processor/references.json` (our
own bibliographic records, from Crossref, PubMed and DataCite metadata).

| Component | What | Licence | Source |
|---|---|---|---|
| Jinja2 3.1.6 (+ MarkupSafe 3.0.4) | The HTML template, with autoescaping. | BSD-3-Clause | https://palletsprojects.com |
| pandas 3.1.0, NumPy 2.5.4 | Tables. | BSD-3-Clause | pandas.pydata.org, numpy.org |
| Matplotlib 3.11.2 | The QC figure. | Matplotlib licence (PSF-based) | https://matplotlib.org |
| contourpy, cycler, fonttools, kiwisolver, packaging, pillow, pyparsing, python-dateutil, six, tzdata | Dependencies of the above. | BSD / MIT / MIT-CMU (Pillow) / Apache-2.0 (tzdata) / Apache-2.0 or BSD-2-Clause (packaging, python-dateutil) | PyPI |
| awslambdaric 4.2.3 (+ simplejson 4.2.0, snapshot-restore-py 1.0.0) | AWS Lambda runtime interface client. | Apache-2.0 (simplejson: MIT / AFL-2.1) | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| python:3.12.16-slim-bookworm | Base image (Debian 12 + CPython 3.12). | PSF (Python); Debian packages under their own free licences | https://hub.docker.com/_/python |

Demo data (not in the image or the repository): Discover dataset 427, Tumanov AV, Shein SA,
Lai Z, Chen Y, Ruparel SB, doi:10.26275/gaap-zlpj, CC BY 4.0 (see `CITATION.cff`).
