# Single-cell pseudobulk project

Replace `data/annotated.h5ad` with an AnnData file whose `.obs` contains the
`sample` and `condition` columns illustrated by `cell-metadata.tsv`. That TSV is
only a guide for preparing `.obs`; the workflow does not consume it directly.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory. Replace the example factor levels
`control` and `treated` in both the AnnData metadata and `workflow.toml`.
