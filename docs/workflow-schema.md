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

- `bulk.rnaseq`, `bulk.salmon`, `bulk.align`, `bulk.rnavar`,
  `bulk.star-reference`, `bulk.de`, and `bulk.enrichment`;
- `single-cell.scrnaseq`, `single-cell.alevin`, `single-cell.scanpy`,
  `single-cell.pseudobulk`, and `single-cell.pseudobulk-de`.

`bulk.salmon` and `single-cell.alevin` are native TxSuite DAGs rather than
nf-core launchers. They build a decoy-aware salmon index or a simpleaf splici
index from `fasta` and `gtf`, or reuse a prebuilt one through `salmon_index`
(with a `tx2gene` table) or `simpleaf_index`. Because TxSuite owns their output
layout, their artifacts resolve at plan time instead of deferring to a
postflight adapter. `bulk.salmon` emits `bulk.gene-counts`, so it substitutes
for `bulk.rnaseq` ahead of `bulk.de`; `single-cell.alevin` emits
`single-cell.matrix`, so it substitutes for `single-cell.scrnaseq` ahead of
`single-cell.scanpy`.

Use `txsuite project validate workflow.toml` instead of relying on this list:
validation is also responsible for artifact types, required parameters,
dependency cycles, backend compatibility, and future registry changes.

## Batch integration

`single-cell.scanpy` can build the neighbor graph, clusters, and UMAP from
Harmony-corrected principal components:

```toml
[workflow.stages.inputs]
input = "${quantify.matrix}"
metadata = "cell-metadata.tsv"   # optional

[workflow.stages.params]
integration = "harmony"
batch_column = "batch"
```

`metadata` is an optional input: supply a cell-metadata TSV keyed by
`barcode_column` when the batch assignment lives outside the matrix, or omit it
when the column is already in the matrix's `obs`. `integration = "harmony"`
requires `batch_column`, and the original PCA, normalized expression, and raw
counts are left untouched.

Harmony is opt-in because it can remove real biology when batch and condition
are confounded. Default is `integration = "none"`.

## Genome alignment

`bulk.align` runs STAR two-pass against a prebuilt index and emits sorted BAMs,
per-sample logs, and a gene-count matrix. STAR produces the counts as a
by-product of alignment (`--quantMode GeneCounts`), so the stage feeds
`bulk.de` directly with no separate quantification step. `star_index` is a
required input, which keeps one implementation of `genomeGenerate` in the repo
and makes reuse of that expensive artifact explicit in the graph.

`ReadsPerGene.out.tab` reports unstranded, forward, and reverse totals side by
side, and the wrong column yields a matrix that is mostly noise without failing.
`strandedness` defaults to `auto`, which classifies the library from the split
between the two *stranded* columns — the unstranded column is approximately
their sum, so comparing all three would always answer `unstranded`. The choice
and its evidence are written to `counts/strandedness.tsv`. The run fails rather
than guessing when the split is neither clearly stranded nor near even, and when
samples disagree with each other. Set `strandedness` explicitly to skip
inference; the report still records what was used.

Samplesheets for `bulk.align` and `bulk.salmon` use the `bulk.fastq-samplesheet`
type: `sample` and `fastq_1` are required, `fastq_2` is optional, and no
`strandedness` column is needed because both stages infer it. That is a weaker
contract than `bulk.rnaseq-samplesheet`, so the two are deliberately distinct
types and cannot be wired interchangeably.

## Variant calling and reference reuse

`bulk.rnavar` launches the pinned nf-core/rnavar release for GATK4 RNA short
variant discovery. Its `variants` output is the `variant_calling/` **directory**
rather than a VCF file, because rnavar writes one VCF per sample and a
file-level contract would match several paths and fail postflight on any
multi-sample run.

`bulk.star-reference` builds a STAR index, FASTA index, and sequence dictionary
once from a genome FASTA and GTF. `bulk.rnavar` accepts `star_index`,
`fasta_fai`, and `dict` as **optional inputs**, so a project can reuse them:

```toml
[[workflow.stages]]
id = "variants"
uses = "bulk.rnavar"
depends_on = ["ref"]

[workflow.stages.inputs]
samplesheet = "samplesheet.csv"
star_index = "${ref.star_index}"
```

Omit them and rnavar derives its own references per run. Optional inputs are
validated exactly like required ones when supplied — artifact types must match
and a reference implies a dependency edge — but a stage plans successfully
without them. Unknown input names are still rejected.

Reference preparation records the read length and splice-junction overhang in
`reference/reference-manifest.tsv`, because an index built for the wrong
overhang runs without complaint while losing junction sensitivity.

The published FASTA index and sequence dictionary keep the caller's FASTA name
— `GRCh38.fa` yields `GRCh38.fa.fai` and `GRCh38.dict` — because GATK resolves
both from the reference basename rather than from the paths it is handed.

TxSuite rejects two rnavar combinations before Nextflow starts: base
recalibration without `dbsnp` or `known_indels`, and an annotation tool without
its cache. Both otherwise surface only after alignment has already run.

## Contrast expansion

`bulk.de` runs one comparison by default. Set `contrasts` to derive the whole
comparison set from the levels of the design column instead of writing one stage
per pair:

```toml
[workflow.stages.params]
design = "condition"
reference = "control"
test = "treated"
contrasts = "vs-reference"    # single (default), vs-reference, or all-pairs
```

`vs-reference` compares every other level against `reference`; `all-pairs`
compares every level combination. `reference` and `test` stay required and name
the *primary* contrast, which keeps its existing outputs: `<method>-results.tsv`,
`significant-genes.tsv`, and the volcano and MA plots are all unchanged. Extra
contrasts are written to `contrasts/DE_<test>_vs_<reference>.tsv`.

Every mode writes a `contrasts.tsv` index — one row in `single` mode — with the
contrast id, its levels, gene and significant-gene counts, and the relative path
to its table. Because the index always exists, the declared output set does not
change with the parameter.

For DESeq2 the expanded contrasts are read off one shared model fit. edgeR and
limma model each comparison as a two-level subset, so those methods refit per
contrast; the tradeoff is that a contrast's numbers are identical whether it was
produced by an expanded run or a single-contrast run.

Above 50 comparisons the stage fails rather than launching the run, naming the
column and its level count. Formula/coefficient mode has no design levels to
expand and rejects any mode other than `single`.

## Artifact references and deferred commands

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
txsuite project init --preset bulk-salmon my-project
txsuite project init --preset bulk-rnavar my-project
txsuite project init --preset scrnaseq my-project
txsuite project init --preset scrnaseq-alevin my-project
txsuite project init --preset scrnaseq-pseudobulk my-project
```

Each preset contains a README, `workflow.toml`, and small syntactically valid
input templates. They are starting points, not runnable biological datasets;
replace every placeholder path and experimental value before validation or run.
