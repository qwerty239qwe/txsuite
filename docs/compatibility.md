# Compatibility

Validated versions for the 0.1 development line:

| Surface | Pinned/tested version | Status |
| --- | --- | --- |
| TxSuite CLI | Python 3.11–3.13 | CI |
| nf-core/rnaseq | 3.26.0 | launcher tested; external raw-data smoke pending |
| nf-core/scrnaseq | 4.2.0 | launcher tested; external raw-data smoke pending |
| Space Ranger | 4.1.0 | external, user-installed; licensed smoke pending |
| Spacemake | 0.9.1b | experimental pass-through |
| DESeq2 | 1.52.0 / Bioconductor 3.23 | Docker smoke passed |
| Scanpy | 1.12.2 / Python 3.12 | Docker smoke passed |
| SpatialData | 0.8.0 | synthetic Docker smoke passed |
| spatialdata-io | 0.7.1 | synthetic Docker smoke passed |
| Squidpy | 1.8.3 | synthetic Docker smoke passed |

“Launcher tested” means validation, command construction, and dry-run behavior
are covered. It does not imply that licensed software or large reference/FASTQ
downloads run in CI.
