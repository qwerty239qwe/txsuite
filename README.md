# TxSuite

Reproducible bulk, single-cell, and spatial transcriptomics from one small CLI.

TxSuite does not reimplement aligners or established pipelines. It provides a
stable configuration, environment checks, and provenance layer over selected
upstream tools.

## Current foundation

```bash
uv run txsuite tools
uv run txsuite tools --modality bulk --stage alignment
uv run txsuite config show
uv run txsuite env doctor
```

For local development:

```bash
python -m pip install -e .
python -m unittest discover -s tests
```

The supported tool matrix and delivery phases are in
[docs/roadmap.md](docs/roadmap.md).

## Project workflows

The project interface plans and runs several typed stages from one versioned
`workflow.toml`. Start from a packaged preset; the target directory must not
already exist:

```bash
txsuite project init --preset bulk-rnaseq my-analysis
cd my-analysis
txsuite project validate workflow.toml
txsuite project plan workflow.toml
txsuite project run workflow.toml --dry-run
txsuite project run workflow.toml
```

`validate` resolves paths, stage contracts, parameters, dependencies, and the
complete command plan without executing tools. `plan` prints the deterministic
plan; add `--json` for machine-readable output. `run` prints the new run-bundle
directory after completion. Inspect it later with:

```bash
txsuite project status results/PROJECT_ID/runs/RUN_ID
```

Select an inclusive range with `--from STAGE --to STAGE`, or a comma-separated
set with `--stages STAGE_A,STAGE_B`. Resume an existing bundle explicitly:

```bash
txsuite project run workflow.toml --resume --run-id RUN_ID
```

Resume skips a completed stage only when its resolved command, parameters,
declared inputs, configuration, and required artifacts still match. Raw
nf-core outputs are resolved after the producer succeeds; downstream symbolic
inputs such as `${rnaseq.counts}` are never passed to a subprocess unresolved.

Available presets are `bulk-rnaseq`, `bulk-salmon`, `scrnaseq`,
`scrnaseq-alevin`, and `scrnaseq-pseudobulk`. The two quantification presets
differ in backend, not in downstream stages: `bulk-rnaseq` and `scrnaseq` launch
pinned nf-core pipelines, while `bulk-salmon` and `scrnaseq-alevin` run native
TxSuite DAGs built on salmon and simpleaf/alevin-fry. The native DAGs quantify
only — no trimming, alignment, or aggregated QC — and their artifact paths are
known at plan time. Replace every placeholder input before running. See the
[workflow schema](docs/workflow-schema.md),
[run-bundle layout](docs/results-layout.md), and
[automation guidance](docs/agent-workflows.md).

## Bulk RNA-seq

TxSuite delegates raw-read QC, trimming, alignment, and quantification to the
pinned nf-core/rnaseq release. Its input CSV requires `sample`, `fastq_1`,
`fastq_2`, and `strandedness` columns.

Templates are provided as `examples/bulk/samplesheet.example.csv` and
`examples/bulk/rnaseq-params.example.json`; replace their placeholder paths.

Preview or launch the workflow:

```bash
uv run txsuite workflow bulk \
  --input samplesheet.csv \
  --outdir results/bulk \
  --params-file rnaseq-params.json \
  --dry-run

uv run txsuite workflow bulk \
  --input samplesheet.csv \
  --outdir results/bulk \
  --params-file rnaseq-params.json
```

Build the owned R image and compare two levels of a metadata column. DESeq2 is
the default; select `--method edger` for edgeR quasi-likelihood or
`--method limma` for limma-voom. Repeat `--covariate` for batch or other
adjustment variables:

```bash
uv run txsuite env build bulk-r

uv run txsuite bulk de \
  --counts counts.tsv \
  --metadata metadata.tsv \
  --design condition \
  --reference control \
  --test treated \
  --covariate batch \
  --padj 0.05 \
  --lfc 1 \
  --outdir results/de

uv run txsuite bulk de \
  --method edger \
  --counts counts.tsv \
  --metadata metadata.tsv \
  --design condition \
  --reference control \
  --test treated \
  --outdir results/edger

# Advanced models use an explicit R formula and one named coefficient.
uv run txsuite bulk de \
  --counts counts.tsv \
  --metadata metadata.tsv \
  --formula "~ batch + condition" \
  --coefficient condition_treated_vs_control \
  --outdir results/adjusted
```

Formula mode is mutually exclusive with `--design`, `--reference`, `--test`,
and `--covariate`. If a coefficient is misspelled, the backend reports the
available names.

All methods write complete and significant DE tables, normalized counts,
library QC, an analysis summary, ordination, MA, volcano, and top-gene plots.
DESeq2 additionally writes VST counts and sample correlations.

Run over-representation analysis or ranked GSEA against any standard GMT gene
set collection:

```bash
uv run txsuite bulk enrich \
  --de results/de/deseq2-results.tsv \
  --genesets pathways.gmt \
  --mode ora \
  --outdir results/ora

uv run txsuite bulk enrich \
  --de results/de/deseq2-results.tsv \
  --genesets pathways.gmt \
  --mode gsea \
  --outdir results/gsea
```

Small synthetic counts, metadata, and GMT inputs are available under
`examples/bulk/` for image smoke tests. Gene identifiers in the DE table and
GMT must use the same namespace.

The counts table must contain genes as rows, samples as columns, and gene IDs in
its first column. Metadata must contain sample IDs in its first column. Every
executed command writes `command.txt`, `stdout.log`, `stderr.log`, `run.json`,
and `txsuite-results.json` under its `.txsuite` run directory.

## Single-cell RNA-seq

Raw reads are handled by pinned nf-core/scrnaseq 4.2.0. Simpleaf is the open
default; `--aligner star` selects STARsolo and `--aligner cellranger` uses a
licensed Cell Ranger installation supplied by the user.

```bash
uv run txsuite workflow single-cell \
  --input samplesheet.csv \
  --outdir results/scrnaseq \
  --params-file scrnaseq-params.json \
  --dry-run

uv run txsuite env build single-cell-python

uv run txsuite single-cell analyze \
  --input examples/single_cell/filtered_feature_bc_matrix \
  --metadata examples/single_cell/metadata.example.tsv \
  --batch-column batch \
  --integration harmony \
  --outdir results/single-cell-smoke \
  --min-genes 1 \
  --min-cells 1 \
  --max-mito-pct 100 \
  --resolution 0.5 \
  --n-hvg 2000 \
  --n-pcs 30 \
  --n-neighbors 15 \
  --doublets score
```

`analyze` accepts a 10x matrix directory, 10x H5 file, or H5AD and writes
`analysis.h5ad`, `cell-qc.tsv`, `clusters.tsv`, `marker-genes.tsv`, a summary,
and the standard TxSuite manifest. Optional metadata is a TSV whose unique
`barcode` column exactly matches the input cells. Scrublet is opt-in with
`--doublets score` or `--doublets filter`; cluster markers use Wilcoxon and are
reported for clusters with at least two cells. They are descriptive rather than
replicate-aware differential-expression results. `--batch-column` makes highly
variable gene selection batch-aware. Add `--integration harmony` to build the
neighbor graph, clusters, and UMAP from Harmony-corrected PCs; the original PCA,
normalized expression, and raw counts remain unchanged. Harmony is opt-in
because it can remove real biology when batch and condition are confounded. A
tiny 10x-format smoke dataset and matching metadata are in `examples/single_cell/`.
Normalization, HVG flavor/count, PCs, neighbors, UMAP distance, marker method,
and the raw-count layer are explicit options. Use `--stop-after qc|pca|clusters`
or `--skip-umap` / `--skip-markers` for partial runs; skipped tabular outputs
remain as header-only files. Existing H5AD input is rebuilt from
`--counts-layer` (default `counts`) so normalized expression is never treated
as raw counts.

For an H5AD whose `obs` contains sample and experimental-design columns, the
native DSL2 workflow runs pseudobulk aggregation and DESeq2 as one resumable
DAG:

```bash
uv run txsuite workflow pseudobulk-de \
  --input annotated.h5ad \
  --sample-column sample \
  --design condition \
  --group-column cell_type \
  --group-value T_cell \
  --covariate batch \
  --reference control \
  --test treated \
  --outdir results/pseudobulk-de \
  --resume
```

The packaged workflow has `local`, `docker`, and `apptainer` profiles and
accepts `--nextflow-config` for cluster settings. Its reusable
`PSEUDOBULK_DE` subworkflow passes comparison metadata through parameter-free
modules; images, resources, failure policy, and publishing remain in
`nextflow.config`. The original
`txsuite single-cell pseudobulk-de` command remains available for running the
same two container steps directly without Nextflow. Group column and value must
be supplied together; one cell type or cluster is analyzed per run. Repeat
`--covariate` for sample-level adjustment variables.

Run many groups or contrasts from a tab-separated manifest with one comparison
per row:

```bash
uv run txsuite single-cell pseudobulk-batch \
  --input annotated.h5ad \
  --sample-column sample \
  --manifest examples/single_cell/comparisons.example.tsv \
  --outdir results/pseudobulk-batch \
  --resume
```

Batch mode uses the same Nextflow DAG by default. Add `--direct` to retain the
sequential Docker runner when Nextflow is unavailable. Required manifest
columns are `comparison`, `design`, `reference`, and `test`.
Optional columns are `group_column`, `group_value`, `method`, comma-separated
`covariates`, `padj`, `lfc`, and `top_genes`. Each comparison keeps its own
outputs under `pseudobulk/<comparison>` and `de/<comparison>`. The workflow
writes `comparison-index.tsv` and `combined-results.tsv` after all runnable
comparisons finish. It continues after individual failures, records failed
comparisons in the index, and exits nonzero when any failed.

### Cell Ranger

`--aligner cellranger` above delegates to nf-core/scrnaseq's own Cell Ranger
module. For a standalone reference build and count, TxSuite also ships a native
resumable `mkref` + `count` Nextflow DAG that calls a user-installed
`cellranger` directly; TxSuite never downloads, bundles, or licenses Cell
Ranger itself.

Build a GRCh38.p14 reference once from GENCODE (verify each SHA-256 against the
[release's published checksum](https://www.gencodegenes.org/human/) before
caching):

```bash
txsuite reference cache \
  --source https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz \
  --sha256 PUBLISHED_64_CHARACTER_SHA256 \
  --name GRCh38.primary_assembly.genome.fa.gz

txsuite reference cache \
  --source https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz \
  --sha256 PUBLISHED_64_CHARACTER_SHA256 \
  --name gencode.v50.primary_assembly.annotation.gtf.gz

gunzip -k .txsuite/references/GRCh38.primary_assembly.genome.fa.gz
gunzip -k .txsuite/references/gencode.v50.primary_assembly.annotation.gtf.gz

uv run txsuite workflow single-cell-cellranger \
  --genome-name GRCh38.p14 \
  --fasta .txsuite/references/GRCh38.primary_assembly.genome.fa \
  --gtf .txsuite/references/gencode.v50.primary_assembly.annotation.gtf \
  --fastqs fastqs/sample1 \
  --sample sample1 \
  --outdir results/cellranger \
  --dry-run
```

Reuse a reference already built by a prior run with `--reference
results/cellranger/reference/GRCh38.p14` instead of `--fasta`/`--gtf`; the DAG
skips `mkref` when a reference path is supplied. `-resume` also skips `mkref`
once its cached output matches. `count` output lands under
`results/cellranger/counts/SAMPLE/outs/filtered_feature_bc_matrix`, in the same
10x layout that `single-cell analyze` accepts.

The `docker` and `apptainer` profiles run this DAG inside a container, but
TxSuite has no image to default to: build one locally from your own
EULA-accepted download and pass its tag through `--cellranger-image`:

```bash
uv run txsuite env build cellranger \
  --tag txsuite/cellranger:local \
  --source-tarball ~/downloads/cellranger-9.0.1.tar.gz
```

`env build cellranger` only packages the install recipe; it copies your
tarball into the build context and never fetches or redistributes it. The
`local` profile instead runs the `cellranger` binary directly from `PATH`.

## Spatial transcriptomics

TxSuite discovers but does not install or redistribute Space Ranger. Preview a
Visium run with:

```bash
uv run txsuite workflow spatial \
  --id sample1 \
  --transcriptome /refs/GRCh38 \
  --fastqs /data/fastqs \
  --sample sample1 \
  --image /data/tissue.tif \
  --slide V19J01-123 \
  --area A1 \
  --outdir results/spaceranger \
  --dry-run
```

Build the owned downstream image, then import a Space Ranger `outs/` directory
into SpatialData and calculate spot QC plus a Squidpy grid-neighbor graph:

```bash
uv run txsuite env build spatial-python
uv run txsuite spatial analyze \
  --input results/spaceranger/sample1/outs \
  --dataset-id sample1 \
  --outdir results/spatial-analysis
```

An existing Spacemake project can be run with the explicitly experimental
`txsuite workflow spatial-open --project-root PROJECT --cores 8` command.

## Reproducibility and HPC

Cache a reference only after verifying its published SHA-256:

```bash
txsuite reference cache \
  --source https://example.org/GRCh38.tar.gz \
  --sha256 PUBLISHED_64_CHARACTER_SHA256 \
  --name GRCh38.tar.gz
```

Release configurations must use immutable `image@sha256:digest` references;
`txsuite env verify-images` reports mutable tags and exits nonzero. See
`docs/slurm.md`, `docs/compatibility.md`, and `docs/release.md` for cluster and
release checks. The `containers` GitHub Actions workflow publishes all four
owned images — `bulk-r`, `salmon`, `single-cell-python`, and `spatial-python` —
to GHCR from version tags or a manual run.

## Project documentation

- [Roadmap and maturity](docs/roadmap.md)
- [Compatibility and current limits](docs/compatibility.md)
- [Workflow schema v1](docs/workflow-schema.md)
- [Run bundles and result provenance](docs/results-layout.md)
- [Agent and automation workflows](docs/agent-workflows.md)
- [Release checklist](docs/release.md)
Owned images intentionally have no fixed `ENTRYPOINT`, so Docker, Nextflow, and
Apptainer can supply an explicit command. Their default `CMD` only prints tool
help, and `/work` is the shared working-directory contract.
