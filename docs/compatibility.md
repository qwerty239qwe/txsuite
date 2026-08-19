# Compatibility

Validated versions for the 0.1 development line:

| Surface | Pinned/tested version | Status |
| --- | --- | --- |
| TxSuite CLI | Python 3.11–3.13 | CI |
| Project workflow schema | v1 | runtime and packaged JSON Schema parity tested |
| Project execution profiles | Docker, Apptainer | schema-supported |
| Project presets | bulk-rnaseq, bulk-salmon, scrnaseq, scrnaseq-alevin, scrnaseq-pseudobulk | packaged-resource and planning tests |
| nf-core/rnaseq | 3.26.0 | launcher tested; external raw-data smoke pending |
| nf-core/scrnaseq | 4.2.0 | launcher tested; external raw-data smoke pending |
| Cell Ranger (native `mkref`+`count` DAG) | user-installed | stub DAG tested in CI; licensed smoke pending |
| Salmon (native `index`+`quant` DAG) | 2.5.1 | stub DAG tested in CI; real FASTQ smoke pending |
| simpleaf (native `index`+`quant` DAG) | 0.28.0 | stub DAG tested in CI; real FASTQ smoke pending |
| alevin-fry | 0.18.0 | USA-mode matrix-market to `.h5ad` conversion smoke passed |
| piscem | 0.22.1 | installed in `txsuite/salmon`; exercised through simpleaf |
| gffread | 0.12.9 | gentrome and transcript-to-gene derivation |
| Space Ranger | 4.1.0 | external, user-installed; licensed smoke pending |
| Spacemake | 0.9.1b | experimental pass-through |
| DESeq2 | 1.52.0 / Bioconductor 3.23 | covariate/formula DE and plot smoke passed |
| edgeR | 4.10.1 / Bioconductor 3.23 | formula-capable quasi-likelihood Docker smoke passed |
| limma | 3.68.4 / Bioconductor 3.23 | formula-capable voom Docker smoke passed |
| clusterProfiler | 4.20.0 / Bioconductor 3.23 | GMT ORA and GSEA smokes passed |
| Scanpy | 1.12.2 / Python 3.12 | count-safe configurable stages, QC, Scrublet, batch-aware HVGs, Leiden markers, and grouped pseudobulk Docker smoke passed |
| Harmonypy | 2.0.0 | PCA integration Docker smoke passed |
| SpatialData | 0.8.0 | synthetic Docker smoke passed |
| spatialdata-io | 0.7.1 | synthetic Docker smoke passed |
| Squidpy | 1.8.3 | synthetic Docker smoke passed |

“Launcher tested” means validation, command construction, and dry-run behavior
are covered. It does not imply that licensed software or large reference/FASTQ
downloads run in CI.

## Project workflow limits

- Schema version 1 is strict: unknown keys, stages, parameters, inputs, outputs,
  dependencies, and artifact-type mismatches are rejected.
- The project schema supports `docker` and `apptainer`. SLURM is configured in
  Nextflow; it is not a third project execution profile.
- `bulk.rnaseq` and `single-cell.scrnaseq` are selected, pinned integrations.
  Their planners, deferred artifacts, and release-specific result adapters are
  tested, but real FASTQ Docker and Apptainer smokes remain release blockers.
- `bulk.de` expands contrasts from the design column when `contrasts` is
  `vs-reference` or `all-pairs`, capped at 50 comparisons. The comparison-set
  logic is unit-tested; the expanded DESeq2, edgeR, and limma runs themselves
  are covered only by the existing single-contrast container smokes.
- `bulk.salmon` and `single-cell.alevin` are native DAGs TxSuite owns end to
  end, so their artifact paths are known at plan time and need no release
  adapter. Their control plane and stubs are tested; real FASTQ runs are not.
- `bulk.de`, `bulk.enrichment`, `single-cell.scanpy`,
  `single-cell.pseudobulk`, and `single-cell.pseudobulk-de` are the currently
  registered downstream stage types. A compatible output type is required for
  every `${stage.artifact}` edge.
- nf-core postflight adapters support the pinned releases listed above. A new
  upstream output layout or release requires a new or updated adapter and an
  explicit compatibility test; TxSuite does not guess across releases.
- Cell Ranger and Space Ranger remain external, user-installed licensed tools.
  Project orchestration does not install them or grant a license. The native
  Cell Ranger `mkref`+`count` DAG (`txsuite workflow single-cell-cellranger`)
  and `env build cellranger` are not part of the schema v1 project workflow
  system, matching how Space Ranger is also invoked outside it.
- Presets contain placeholders and tiny syntax fixtures, not reference genomes,
  FASTQs, or validated biological examples.

## Platform and packaging limits

The dependency-free control plane is tested on CPython 3.11–3.13 in GitHub
Actions. Analysis backends have their own Linux/container requirements. The
wheel contains the workflow JSON Schema and all preset resources and is audited
in CI. Windows can load, validate, plan, and scaffold projects, but the selected
container, Nextflow, and licensed-tool execution paths are not claimed as
end-to-end validated there.
