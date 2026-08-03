# Single-cell RNA-seq project

Replace the `data/...` FASTQ placeholders and review the genome and protocol in
`scrnaseq-params.json` and `workflow.toml`.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory. The Scanpy thresholds are valid starting
values, not study-specific recommendations.
