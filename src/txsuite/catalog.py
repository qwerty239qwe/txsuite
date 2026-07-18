from __future__ import annotations

from typing import NamedTuple


class Tool(NamedTuple):
    modality: str
    stage: str
    name: str
    executable: str
    environment: str
    status: str


TOOLS = (
    Tool(
        "bulk",
        "workflow",
        "nf-core/rnaseq",
        "nextflow",
        "upstream containers",
        "selected",
    ),
    Tool("bulk", "qc", "FastQC", "fastqc", "nf-core/rnaseq", "selected"),
    Tool("bulk", "qc", "MultiQC", "multiqc", "nf-core/rnaseq", "selected"),
    Tool("bulk", "trim", "fastp", "fastp", "nf-core/rnaseq", "selected"),
    Tool("bulk", "alignment", "STAR", "STAR", "nf-core/rnaseq", "selected"),
    Tool("bulk", "quantification", "Salmon", "salmon", "nf-core/rnaseq", "selected"),
    Tool(
        "bulk",
        "differential-expression",
        "DESeq2",
        "docker",
        "txsuite/bulk-r",
        "ready",
    ),
    Tool(
        "bulk",
        "enrichment",
        "clusterProfiler ORA/GSEA",
        "docker",
        "txsuite/bulk-r",
        "ready",
    ),
    Tool(
        "single-cell",
        "workflow",
        "nf-core/scrnaseq",
        "nextflow",
        "upstream containers",
        "selected",
    ),
    Tool("single-cell", "qc", "FastQC", "fastqc", "nf-core/scrnaseq", "selected"),
    Tool(
        "single-cell",
        "qc",
        "Scanpy",
        "docker",
        "txsuite/single-cell-python",
        "ready",
    ),
    Tool(
        "single-cell",
        "alignment",
        "Simpleaf",
        "simpleaf",
        "nf-core/scrnaseq",
        "selected",
    ),
    Tool(
        "single-cell", "alignment", "STARsolo", "STAR", "nf-core/scrnaseq", "selected"
    ),
    Tool(
        "single-cell",
        "workflow",
        "Cell Ranger",
        "cellranger",
        "user install",
        "external",
    ),
    Tool(
        "single-cell",
        "differential-expression",
        "pseudobulk DESeq2",
        "docker",
        "txsuite/single-cell-python + bulk-r",
        "ready",
    ),
    Tool(
        "spatial", "workflow", "Space Ranger", "spaceranger", "user install", "external"
    ),
    Tool(
        "spatial",
        "qc",
        "SpatialData + Squidpy",
        "docker",
        "txsuite/spatial-python",
        "ready",
    ),
    Tool("spatial", "workflow", "Spacemake", "spacemake", "Apptainer", "experimental"),
)


def select_tools(
    modality: str | None = None, stage: str | None = None
) -> tuple[Tool, ...]:
    return tuple(
        tool
        for tool in TOOLS
        if (modality is None or tool.modality == modality)
        and (stage is None or tool.stage == stage)
    )
