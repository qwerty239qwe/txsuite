from __future__ import annotations

import csv
import re
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError, run_command


REQUIRED_COLUMNS = {"sample", "fastq_1", "fastq_2", "strandedness"}
STRANDEDNESS = {"auto", "forward", "reverse", "unstranded"}


def validate_samplesheet(path: Path) -> int:
    if not path.is_file():
        raise TxSuiteError(f"Samplesheet does not exist: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or ())
        if missing:
            raise TxSuiteError(
                f"Samplesheet is missing columns: {', '.join(sorted(missing))}"
            )
        rows = 0
        for line, row in enumerate(reader, start=2):
            if (
                not (row.get("sample") or "").strip()
                or not (row.get("fastq_1") or "").strip()
            ):
                raise TxSuiteError(
                    f"Samplesheet line {line} requires sample and fastq_1"
                )
            strandedness = (row.get("strandedness") or "").strip().lower()
            if strandedness not in STRANDEDNESS:
                raise TxSuiteError(
                    f"Samplesheet line {line} has invalid strandedness: {row['strandedness']}"
                )
            rows += 1
    if not rows:
        raise TxSuiteError("Samplesheet has no data rows")
    return rows


def workflow_command(
    config: dict[str, Any],
    *,
    samplesheet: Path,
    outdir: Path,
    params_file: Path | None = None,
    nextflow_config: Path | None = None,
    resume: bool = False,
) -> list[str]:
    pipeline = config["pipelines"]["bulk"]
    command = [
        "nextflow",
        "run",
        pipeline["name"],
        "-r",
        pipeline["release"],
        "-profile",
        config["execution"]["profile"],
        "--input",
        str(samplesheet.resolve()),
        "--outdir",
        str(outdir.resolve()),
    ]
    if params_file is not None:
        if not params_file.is_file():
            raise TxSuiteError(f"Params file does not exist: {params_file}")
        command.extend(["-params-file", str(params_file.resolve())])
    if nextflow_config is not None:
        if not nextflow_config.is_file():
            raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")
        command.extend(["-c", str(nextflow_config.resolve())])
    if resume:
        command.append("-resume")
    return command


def deseq2_command(
    *,
    image: str,
    counts: Path,
    metadata: Path,
    design: str,
    reference: str,
    test: str,
    outdir: Path,
    check_inputs: bool = True,
) -> list[str]:
    if check_inputs:
        for label, path in (("Counts", counts), ("Metadata", metadata)):
            if not path.is_file():
                raise TxSuiteError(f"{label} file does not exist: {path}")
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", design):
        raise TxSuiteError("Design must be a simple metadata column name")
    if not reference or not test or reference == test:
        raise TxSuiteError("Reference and test levels must be non-empty and different")
    if not image.strip():
        raise TxSuiteError("DESeq2 image cannot be empty")
    mounts = (
        f"type=bind,source={counts.resolve()},target=/input/counts.tsv,readonly",
        f"type=bind,source={metadata.resolve()},target=/input/metadata.tsv,readonly",
        f"type=bind,source={outdir.resolve()},target=/output",
    )
    return [
        "docker",
        "run",
        "--rm",
        "--mount",
        mounts[0],
        "--mount",
        mounts[1],
        "--mount",
        mounts[2],
        image,
        "/input/counts.tsv",
        "/input/metadata.tsv",
        design,
        reference,
        test,
        "/output",
    ]


def build_bulk_r_image(tag: str, *, run_dir: Path) -> None:
    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    package = resources.files("txsuite.resources.bulk_r")
    with tempfile.TemporaryDirectory(prefix="txsuite-bulk-r-") as directory:
        context = Path(directory)
        for name in ("Dockerfile", "deseq2.R"):
            (context / name).write_text(
                package.joinpath(name).read_text(encoding="utf-8"), encoding="utf-8"
            )
        run_command(
            ["docker", "build", "--tag", tag, str(context)],
            run_dir=run_dir,
            task="env.build.bulk-r",
            backend="docker",
            inputs={
                "base": (
                    "bioconductor/bioconductor_docker:RELEASE_3_23@"
                    "sha256:1d871e1ca9cca76b220eb16e22677e728f4352f81a9ee91aaf29e24aea43e624"
                )
            },
            outputs={"image": tag},
            artifacts=[{"kind": "container-image", "label": "bulk-r", "path": tag}],
        )
