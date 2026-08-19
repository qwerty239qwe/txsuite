"""Splice-aware genome alignment for bulk RNA-seq.

STAR produces the gene counts as a by-product of alignment
(``--quantMode GeneCounts``), so this stage yields both sorted alignments and a
count matrix that drops straight into ``bulk.de`` without a second tool.

The STAR index is a required input rather than something rebuilt here: it is
the expensive artifact ``bulk.star-reference`` exists to produce, and keeping
one implementation of ``genomeGenerate`` in the repo keeps its parameters --
notably the splice-junction overhang -- in a single place.
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Any

from txsuite.bulk import validate_quant_samplesheet
from txsuite.runtime import TxSuiteError

STRANDEDNESS = ("auto", "unstranded", "forward", "reverse")


def align_workflow_command(
    config: dict[str, Any],
    *,
    samplesheet: Path,
    star_index: Path,
    outdir: Path,
    strandedness: str = "auto",
    two_pass: bool = True,
    star_image: str | None = None,
    threads: int = 4,
    memory_gb: int = 32,
    nextflow_config: Path | None = None,
    resume: bool = False,
    check_inputs: bool = True,
) -> list[str]:
    """Build the native STAR alignment and gene-counting Nextflow DAG."""

    if strandedness not in STRANDEDNESS:
        raise TxSuiteError(
            f"Strandedness must be one of: {', '.join(STRANDEDNESS)}"
        )
    if threads < 1 or memory_gb < 1:
        raise TxSuiteError("STAR threads and memory must be positive")
    if check_inputs:
        validate_quant_samplesheet(samplesheet)
        if not star_index.is_dir():
            raise TxSuiteError(f"STAR index does not exist: {star_index}")
    if nextflow_config is not None and not nextflow_config.is_file():
        raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")

    image = star_image or config["images"]["star"]
    if not image.strip():
        raise TxSuiteError("Workflow images cannot be empty")

    workflow = resources.files("txsuite.resources.nextflow").joinpath("bulk_align.nf")
    command = [
        "nextflow",
        "run",
        str(workflow),
        "-profile",
        config["execution"]["profile"],
        "-work-dir",
        str((outdir / ".nextflow-work").resolve()),
    ]
    if nextflow_config is not None:
        command.extend(["-c", str(nextflow_config.resolve())])
    if resume:
        command.append("-resume")
    command.extend(
        [
            "--samplesheet",
            str(samplesheet.resolve()),
            "--star_index",
            str(star_index.resolve()),
            "--outdir",
            str(outdir.resolve()),
            "--star_image",
            image,
            "--star_strandedness",
            strandedness,
            "--star_two_pass",
            "true" if two_pass else "false",
            "--star_threads",
            str(threads),
            "--star_memory_gb",
            str(memory_gb),
        ]
    )
    return command


__all__ = ["STRANDEDNESS", "align_workflow_command"]
