# SLURM

Nextflow submits individual pipeline processes to SLURM; do not wrap the whole
pipeline in a custom TxSuite scheduler. Start from
`examples/slurm/nextflow.config`, then change its queue and resource limits to
match the cluster.

Set the TxSuite execution profile to Apptainer:

```toml
[execution]
profile = "apptainer"
```

Then launch either nf-core workflow with the cluster configuration:

```bash
txsuite workflow bulk \
  --input samplesheet.csv \
  --outdir results/bulk \
  --params-file rnaseq-params.json \
  --nextflow-config examples/slurm/nextflow.config \
  --resume
```

Use an institutional profile from `nf-core/configs` when one exists. The example
is intentionally local: account names, partitions, bind paths, and limits are
site policy and do not belong in TxSuite defaults.

Space Ranger and Spacemake are external applications. Run them according to the
cluster administrator's supported module/container and submission policy.
