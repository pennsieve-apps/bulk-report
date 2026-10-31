# sample-input

bulk-report reads bulk-de's output (port `de`, one `<LABEL>/` folder per comparison) and,
optionally, rnaseq-quant's counts output (port `counts`) for read-level QC.

## Synthetic fixture (no download)

```sh
make fixture PYTHON=.venv/bin/python   # data/fixture/{linked,flat}/
make run PYTHON=.venv/bin/python       # Docker over data/fixture/linked -> data/output/index.html
```

`linked/de/condition/` and `linked/de/sex-only/` are two bulk-de runs, `linked/counts/` an
rnaseq-quant output (one library with a 31 % mapping rate, flagged); `flat/` holds the same
without port folders. The library names are Discover 427's, so its data credit is recognised.

## Real data

```sh
docker run --rm --platform linux/amd64 --memory 2g --network none --read-only --tmpfs /tmp \
  --env-file dev.env -v /path/to/bulk-de/output:/data/input/de:ro \
  -v /path/to/rnaseq-quant/counts:/data/input/counts:ro -v "$PWD/data/output:/data/output" bulk-report:dev
```
