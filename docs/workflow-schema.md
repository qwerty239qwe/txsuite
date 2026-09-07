# Project workflow schema v1

TxSuite project files are TOML documents validated against schema version 1.
The packaged JSON Schema is
`txsuite.resources.project/workflow.schema.json`; it is also included in the
wheel for editors and external validators. TxSuite's runtime validation is the
authoritative interpretation of TOML types and stage registry contracts.

## Top-level contract

```toml
schema_version = 1

[project]
id = "analysis_id"
modality = "bulk"                 # bulk, single-cell, or spatial
output_root = "results"

[execution]
profile = "docker"                # docker or apptainer
resume = true

[[workflow.stages]]
id = "rnaseq"
uses = "bulk.rnaseq"

[workflow.stages.inputs]
samplesheet = "samplesheet.csv"
```

Unknown keys are errors at every schema layer. Project and stage IDs begin with
a letter and contain only letters, numbers, `_`, or `-`. Relative paths resolve
against the directory containing `workflow.toml`, except tool-specific parameter
files whose command builders document a working-directory convention.

`project.modality` must match every selected stage. Version 1 execution profiles
are `docker` and `apptainer`; scheduler integration is supplied through a
Nextflow configuration rather than another project profile.

## Stages

Each `[[workflow.stages]]` has:

- required `id` and registered `uses` values;
- optional unique `depends_on` stage IDs;
- typed `inputs` matching the registry exactly;
- stage-specific `params`, with defaults applied by the registry;
- optional named `outputs` that override planned artifact paths.

Currently registered stage types are:

- `bulk.rnaseq`, `bulk.de`, and `bulk.enrichment`;
- `single-cell.scrnaseq`, `single-cell.scanpy`,
  `single-cell.pseudobulk`, and `single-cell.pseudobulk-de`.

Use `txsuite project validate workflow.toml` instead of relying on this list:
validation is also responsible for artifact types, required parameters,
dependency cycles, backend compatibility, and future registry changes.

## Artifact references and deferred commands

### Single-cell controls

`single-cell.scanpy` accepts the analysis CLI controls as snake_case parameters:
`metadata`, `barcode_column`, `batch_column`, `integration`, `counts_layer`,
`target_sum`, `n_hvg`, `hvg_flavor`, `n_pcs`, `n_neighbors`, `umap_min_dist`,
`marker_method`, `stop_after`, `skip_umap`, `skip_markers`, `doublets`,
`doublet_batch_column`, `expected_doublet_rate`, `doublet_threshold`, and
`top_markers`, in addition to its QC thresholds and resolution.

Both pseudobulk stages accept `counts_layer`, `group_column`, `group_value`,
and `covariates` (an array of column names). Group column and value must be
provided together. The native DE stage publishes its `de_results` artifact at
`de/<test>_vs_<reference>/deseq2-results.tsv` beneath the stage output directory.

Metadata and configuration file parameters are resolved relative to the
workflow file and fingerprinted for resume. Change their contents to invalidate
the affected stage even when the filename stays the same.

### Referencing outputs

An input can reference a producer's named output using the exact form
`${stage_id.artifact_id}`:

```toml
[[workflow.stages]]
id = "differential"
uses = "bulk.de"
depends_on = ["rnaseq"]

[workflow.stages.inputs]
counts = "${rnaseq.counts}"
metadata = "metadata.tsv"
```

The referenced producer must exist, expose that artifact name, and produce the
type required by the consumer. A reference also implies a dependency even when
it is omitted from `depends_on`.

Some nf-core paths are unknowable until the raw workflow finishes. Their
commands are marked deferred in the plan. After the producer exits successfully,
TxSuite applies the adapter for the pinned pipeline release and requires one
matching path for every declared dynamic artifact. Zero or multiple matches fail
the producer unless an explicit output override selects an existing path. Only
then does TxSuite substitute the absolute path into each downstream argv element.
It never invokes a command containing an unresolved project artifact token.

## Presets

Create a safe, non-overwriting scaffold with:

```bash
txsuite project init --preset bulk-rnaseq my-project
txsuite project init --preset scrnaseq my-project
txsuite project init --preset scrnaseq-pseudobulk my-project
```

Each preset contains a README, `workflow.toml`, and small syntactically valid
input templates. They are starting points, not runnable biological datasets;
replace every placeholder path and experimental value before validation or run.
