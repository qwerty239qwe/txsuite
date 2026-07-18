# TxSuite roadmap

## Goal

One configuration should reproduce the same transcriptomics analysis on a
workstation or HPC system. TxSuite owns configuration, environment selection,
command validation, provenance, and stable result manifests. Upstream tools own
the biological algorithms.

## Non-goals

- Reimplementing QC, alignment, quantification, clustering, or statistics.
- Maintaining forks of nf-core pipelines.
- Downloading or redistributing licensed 10x Genomics software.
- Supporting every assay and backend in the first release.

## Product surface

1. **Workflow API**: reproducible raw-data pipelines.
2. **CLI API**: individual QC, alignment, and analysis tasks.
3. **Python API**: downstream operations on AnnData and SpatialData objects.

All three APIs will use the same configuration and emit the same run manifest.

## Environment ownership

- nf-core workflows use their own pinned process containers; TxSuite does not
  rebuild them.
- TxSuite builds only the downstream images it owns: `bulk-r`,
  `single-cell-python`, and `spatial-python`.
- Cell Ranger and Space Ranger are installed and licensed by the user; TxSuite
  only discovers and invokes them.
- `txsuite env doctor`, `env list`, and `env build` form the environment API.
- Docker is the workstation default and Apptainer is the HPC target.

## Configuration contract

Configuration precedence is:

1. packaged defaults;
2. user config (`~/.config/txsuite/config.toml`, or `%APPDATA%/txsuite/config.toml`);
3. project `txsuite.toml`;
4. explicit CLI options.

Pipeline releases and container image digests must be recorded in every run.
Secrets and machine-specific paths do not belong in committed project config.

## Selected raw-data tools

| Modality | Stage | Default | Alternatives | Delivery |
| --- | --- | --- | --- | --- |
| Bulk RNA-seq | QC | FastQC + MultiQC | fastp reports | nf-core/rnaseq |
| Bulk RNA-seq | Alignment | STAR | HISAT2 | nf-core/rnaseq |
| Bulk RNA-seq | Quantification | Salmon | RSEM, featureCounts | nf-core/rnaseq |
| Single-cell | QC | FastQC + pipeline metrics | CellBender later | nf-core/scrnaseq |
| Single-cell | Alignment/counting | Simpleaf | STARsolo, Kallisto/BUS | nf-core/scrnaseq |
| Single-cell 10x | Alignment/counting | user-installed Cell Ranger | Simpleaf | external licensed backend |
| Spatial Visium | QC/alignment/counting | user-installed Space Ranger | Spacemake later | external licensed backend |
| Other sequencing-based spatial | Processing | Spacemake | assay-specific pipeline | experimental, phase 3 |

## Delivery phases

### Phase 0 — foundation

Deliver:

- dependency-free Python CLI;
- layered TOML configuration;
- curated tool catalog;
- environment diagnostics;
- CI running the stdlib test suite.

Done when `txsuite tools`, `txsuite config show`, and `txsuite env doctor` work
without installing analysis libraries.

### Phase 1 — bulk RNA-seq

Implementation status: the launcher, validation, and provenance layer are
implemented. The pinned DESeq2 1.52.0 image passed a Docker smoke analysis.
Docker and Apptainer nf-core end-to-end dataset validation is still required
before the raw workflow backend is labelled **ready**.

Deliver:

- `txsuite workflow bulk` launching pinned `nf-core/rnaseq`;
- sample-sheet validation;
- Docker and Apptainer profiles;
- `txsuite env build bulk-r` for the owned DESeq2 image;
- provenance and result manifest;
- task commands for FastQC, fastp, STAR, and Salmon only when they provide value
  outside the full workflow;
- downstream DESeq2 through a separate R image.

Done when a small paired-end test dataset produces counts, MultiQC output, and a
complete manifest under both Docker and Apptainer.

### Phase 2 — single-cell RNA-seq

Implementation status: the nf-core/scrnaseq launcher, Scanpy image, standard
analysis path, and pseudobulk-to-DESeq2 path are implemented. The owned image
has passed the bundled 10x smoke dataset; raw FASTQ end-to-end validation remains
an external integration test.

Deliver:

- `txsuite workflow single-cell` launching pinned `nf-core/scrnaseq`;
- Simpleaf as the open default and STARsolo as an alternative;
- optional discovery of a user-installed Cell Ranger;
- `txsuite env build single-cell-python` for the owned Scanpy image;
- AnnData import plus Scanpy QC, normalization, PCA, neighbors, Leiden, and UMAP;
- pseudobulk differential expression through the bulk DESeq2 backend.

Done when a small 10x dataset produces H5AD, QC metrics, clusters, and a manifest.

### Phase 3 — spatial transcriptomics

Implementation status: the Space Ranger launcher, experimental Spacemake
pass-through, and owned SpatialData/Squidpy analysis image are implemented.
The owned image passed a synthetic SpatialData smoke test. A public Visium
download and a licensed Space Ranger end-to-end run remain external checks.

Deliver:

- Visium validation and a user-installed Space Ranger launcher;
- `txsuite env build spatial-python` for the owned SpatialData/Squidpy image;
- SpatialData import and basic Squidpy QC/neighbor analysis;
- Spacemake as an explicitly experimental open backend;
- no Xenium/CosMx/MERSCOPE raw processing until a stable upstream workflow is
  selected.

Done when a public Visium sample produces a SpatialData object, QC report, and
manifest without TxSuite redistributing licensed software.

### Phase 4 — hardening

Implementation status: checksum-verified reference caching, immutable image
validation, native resume guidance, a custom Nextflow config hook, SLURM example,
release/compatibility documentation, and GHCR publishing automation are
implemented. Published image digests and the external raw-data smoke matrix
remain release blockers.

Deliver only after the three modality smoke tests pass:

- resumable runs;
- reference cache with checksums;
- image digest locking;
- SLURM examples;
- release documentation and compatibility table.

## Backend maturity labels

- **ready**: wrapped, pinned, tested, and documented.
- **selected**: chosen for implementation but not yet runnable through TxSuite.
- **experimental**: runnable only with explicit opt-in and version warnings.
- **external**: discovered and invoked, but installed/licensed by the user.

Nothing is advertised as ready until a small end-to-end dataset passes in CI or
documented external integration testing.

## Upstream references

- [nf-core/rnaseq](https://nf-co.re/rnaseq/latest/)
- [nf-core/scrnaseq](https://nf-co.re/scrnaseq/latest/)
- [10x Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/latest)
- [Spacemake](https://spacemake.readthedocs.io/en/latest/)
