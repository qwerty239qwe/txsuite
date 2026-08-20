"""Reference index preparation for splice-aware alignment and GATK.

This module builds indexes from a genome FASTA and annotation GTF.  It is the
counterpart to :mod:`txsuite.hardening`, which acquires and checksum-verifies
reference *downloads*; nothing here fetches anything from the network.

The artifacts are deliberately shared: a STAR index serves both alignment and
nf-core/rnavar, and the FASTA index and sequence dictionary that GATK requires
are derived from the same FASTA in the same run.
"""

from __future__ import annotations

import re
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError, run_command


def star_reference_workflow_command(
    config: dict[str, Any],
    *,
    fasta: Path,
    gtf: Path,
    outdir: Path,
    read_length: int = 100,
    sjdb_overhang: int | None = None,
    star_image: str | None = None,
    threads: int = 4,
    memory_gb: int = 32,
    genome_sa_index_nbases: int | None = None,
    nextflow_config: Path | None = None,
    resume: bool = False,
    check_inputs: bool = True,
) -> list[str]:
    """Build the native STAR/GATK reference-preparation Nextflow DAG.

    ``sjdb_overhang`` defaults to ``read_length - 1``, which is what STAR
    documents for splice-junction database construction. An index built for the
    wrong overhang still runs but loses junction sensitivity, so the effective
    value is recorded in the published manifest rather than left implicit.
    """

    if read_length < 2:
        raise TxSuiteError("Read length must be at least 2 to derive an overhang")
    overhang = read_length - 1 if sjdb_overhang is None else sjdb_overhang
    if overhang < 1:
        raise TxSuiteError("Splice-junction overhang must be positive")
    if threads < 1 or memory_gb < 1:
        raise TxSuiteError("STAR threads and memory must be positive")
    if genome_sa_index_nbases is not None and not 1 <= genome_sa_index_nbases <= 16:
        raise TxSuiteError("STAR genomeSAindexNbases must be between 1 and 16")
    if check_inputs:
        if not fasta.is_file():
            raise TxSuiteError(f"Genome FASTA does not exist: {fasta}")
        if not gtf.is_file():
            raise TxSuiteError(f"Annotation GTF does not exist: {gtf}")
    if nextflow_config is not None and not nextflow_config.is_file():
        raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")

    image = star_image or config["images"]["star"]
    if not image.strip():
        raise TxSuiteError("Workflow images cannot be empty")

    workflow = resources.files("txsuite.resources.nextflow").joinpath(
        "bulk_star_reference.nf"
    )
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
            "--fasta",
            str(fasta.resolve()),
            "--gtf",
            str(gtf.resolve()),
            "--outdir",
            str(outdir.resolve()),
            "--star_image",
            image,
            "--star_read_length",
            str(read_length),
            "--star_sjdb_overhang",
            str(overhang),
            "--star_threads",
            str(threads),
            "--star_memory_gb",
            str(memory_gb),
        ]
    )
    if genome_sa_index_nbases is not None:
        command.extend(["--star_sa_index_nbases", str(genome_sa_index_nbases)])
    return command


def build_star_image(tag: str, *, run_dir: Path) -> None:
    """Build the reference-preparation image from the packaged recipe."""

    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._:/-]*", tag):
        raise TxSuiteError("Image tag contains unsupported characters")
    package = resources.files("txsuite.resources.star")
    with tempfile.TemporaryDirectory(prefix="txsuite-star-") as directory:
        context = Path(directory)
        (context / "Dockerfile").write_text(
            package.joinpath("Dockerfile").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        run_command(
            ["docker", "build", "--tag", tag, str(context)],
            run_dir=run_dir,
            task="env.build.star",
            backend="docker",
            inputs={
                "base": (
                    "condaforge/miniforge3:26.3.2-3@"
                    "sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9"
                ),
                "STAR": "2.7.11b",
                "samtools": "1.24",
                "gatk4": "4.6.2.0",
            },
            outputs={"image": tag},
            artifacts=[{"kind": "container-image", "label": "star", "path": tag}],
        )


__all__ = ["build_star_image", "star_reference_workflow_command"]
