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
```

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
  --outdir results/single-cell-smoke \
  --min-genes 1 \
  --min-cells 1 \
  --max-mito-pct 100 \
  --resolution 0.5
```

`analyze` accepts a 10x matrix directory, 10x H5 file, or H5AD and writes
`analysis.h5ad`, `cell-qc.tsv`, `clusters.tsv`, a summary, and the standard
TxSuite manifest. A tiny 10x-format smoke dataset is in `examples/single_cell/`.

For an H5AD whose `obs` contains sample and experimental-design columns, the
native DSL2 workflow runs pseudobulk aggregation and DESeq2 as one resumable
DAG:

```bash
uv run txsuite workflow pseudobulk-de \
  --input annotated.h5ad \
  --sample-column sample \
  --design condition \
  --reference control \
  --test treated \
  --outdir results/pseudobulk-de \
  --resume
```

The packaged workflow has `local`, `docker`, and `apptainer` profiles and
accepts `--nextflow-config` for cluster settings. The original
`txsuite single-cell pseudobulk-de` command remains available for running the
same two container steps directly without Nextflow.

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
release checks. The `containers` GitHub Actions workflow publishes all three
owned images to GHCR from version tags or a manual run.
