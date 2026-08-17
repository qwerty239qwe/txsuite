# Single-cell alevin project

Droplet quantification with simpleaf and alevin-fry, without an nf-core
pipeline. The DAG builds a splici index from `data/genome.fa` and
`data/genes.gtf`, quantifies every sample, and writes one `.h5ad` whose
`sample` column labels each barcode, ready for `single-cell.scanpy`.

`fastq_1` must hold the barcode/UMI read and `fastq_2` the cDNA read; that is
the standard 10x layout. Replace every `data/...` path with real inputs.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory. Set `chemistry` to match the kit that
generated the libraries, and pass `whitelist` when you have a permit list of
your own; without it the run uses the unfiltered barcode list. To reuse an
index built earlier, drop `fasta` and `gtf` and set `simpleaf_index` instead.

The Scanpy thresholds are valid starting values, not study-specific
recommendations.
