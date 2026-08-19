# Bulk salmon project

Selective-alignment quantification without an aligner or an nf-core pipeline.
The DAG builds a decoy-aware index from `data/genome.fa` and `data/genes.gtf`,
quantifies every sample, and summarizes transcripts to a gene-count matrix that
feeds `bulk.de` directly.

Replace every `data/...` path with real inputs. Keep `samplesheet.csv` sample
IDs identical to the `sample` column in `metadata.tsv`.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory. To reuse an index built earlier, drop
`fasta` and `gtf` from the `quantify` params and set both `salmon_index` and
`tx2gene` instead; a prebuilt index carries no annotation, so the table is
required for gene-level counts.

The values `control` and `treated` are valid placeholders; change them together
in `metadata.tsv` and `workflow.toml`.
