# Bulk RNA-seq variant project

GATK short variant discovery from RNA-seq through the pinned nf-core/rnavar
release: STAR two-pass alignment, MarkDuplicates, SplitNCigarReads, base
recalibration, HaplotypeCaller, and hard filtering.

The `ref` stage builds the STAR index, FASTA index, and sequence dictionary
once, and `variants` reuses them through `${ref.star_index}`. Those inputs are
optional: delete them and rnavar derives its own references per run, which is
simpler but rebuilds a STAR index every time — around an hour and 32 GB of RAM
for a human genome.

Replace every `data/...` path with real inputs.

```console
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml
```

Run these commands from this directory.

Two settings deserve attention before a real run:

- `read_length` drives STAR's splice-junction overhang (`read_length - 1`). Set
  it to your actual read length; an index built for the wrong value still runs
  but loses junction sensitivity. The value used is recorded in
  `reference/reference-manifest.tsv`.
- `dbsnp` supplies known sites for base recalibration. Without known sites, set
  `skip_baserecalibration = true` instead — TxSuite refuses the combination of
  recalibration and no known sites rather than failing hours into the run.

Variants land in `variant_calling/<sample>/`, one VCF per sample. To annotate
them, set `tools = ["snpeff"]` or `["vep"]` and supply the matching
`snpeff_cache` or `vep_cache` directory.
