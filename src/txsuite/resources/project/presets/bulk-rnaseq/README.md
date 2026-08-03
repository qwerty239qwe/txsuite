# Bulk RNA-seq project

Replace every `data/...` path and the example gene set with real inputs. Keep
`samplesheet.csv` sample IDs identical to the `sample` column in `metadata.tsv`.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory. The values `control` and `treated` are
valid placeholders; change them together in `metadata.tsv` and `workflow.toml`.
