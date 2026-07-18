from __future__ import annotations

import csv
import re
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError, run_command


REQUIRED_COLUMNS = ("sample", "fastq_1", "fastq_2")
ALIGNERS = ("simpleaf", "star", "cellranger")


def validate_samplesheet(path: Path) -> int:
    if not path.is_file():
        raise TxSuiteError(f"Samplesheet does not exist: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple((reader.fieldnames or ())[:3]) != REQUIRED_COLUMNS:
            raise TxSuiteError(
                "Samplesheet's first columns must be sample, fastq_1, fastq_2"
            )
        rows = 0
        for line, row in enumerate(reader, start=2):
            if any(not (row.get(column) or "").strip() for column in REQUIRED_COLUMNS):
                raise TxSuiteError(
                    f"Samplesheet line {line} requires sample, fastq_1, and fastq_2"
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
    aligner: str = "simpleaf",
    protocol: str | None = None,
    params_file: Path | None = None,
    nextflow_config: Path | None = None,
    resume: bool = False,
) -> list[str]:
    if aligner not in ALIGNERS:
        raise TxSuiteError(f"Aligner must be one of: {', '.join(ALIGNERS)}")
    pipeline = config["pipelines"]["single_cell"]
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
        "--aligner",
        aligner,
    ]
    if protocol:
        command.extend(["--protocol", protocol])
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


def analysis_command(
    *,
    image: str,
    input_path: Path,
    outdir: Path,
    min_genes: int,
    min_cells: int,
    max_mito_pct: float,
    resolution: float,
) -> list[str]:
    if not input_path.exists():
        raise TxSuiteError(f"Single-cell input does not exist: {input_path}")
    if not image.strip():
        raise TxSuiteError("Single-cell image cannot be empty")
    if min_genes < 0 or min_cells < 0 or not 0 <= max_mito_pct <= 100:
        raise TxSuiteError(
            "QC thresholds must be non-negative; max mito percent <= 100"
        )
    if resolution <= 0:
        raise TxSuiteError("Leiden resolution must be positive")
    target = "/input/data" if input_path.is_dir() else f"/input/{input_path.name}"
    return [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={input_path.resolve()},target={target},readonly",
        "--mount",
        f"type=bind,source={outdir.resolve()},target=/output",
        image,
        "analyze",
        target,
        "/output",
        "--min-genes",
        str(min_genes),
        "--min-cells",
        str(min_cells),
        "--max-mito-pct",
        str(max_mito_pct),
        "--resolution",
        str(resolution),
    ]


def pseudobulk_command(
    *, image: str, h5ad: Path, outdir: Path, sample_column: str, design: str
) -> list[str]:
    if not h5ad.is_file():
        raise TxSuiteError(f"H5AD file does not exist: {h5ad}")
    if not image.strip():
        raise TxSuiteError("Single-cell image cannot be empty")
    for label, value in (("Sample column", sample_column), ("Design", design)):
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", value):
            raise TxSuiteError(f"{label} must be a simple column name")
    return [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={h5ad.resolve()},target=/input/data.h5ad,readonly",
        "--mount",
        f"type=bind,source={outdir.resolve()},target=/output",
        image,
        "pseudobulk",
        "/input/data.h5ad",
        "/output",
        "--sample-column",
        sample_column,
        "--design",
        design,
    ]


def build_single_cell_image(tag: str, *, run_dir: Path) -> None:
    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    package = resources.files("txsuite.resources.single_cell_python")
    with tempfile.TemporaryDirectory(prefix="txsuite-single-cell-") as directory:
        context = Path(directory)
        for name in ("Dockerfile", "single_cell.py"):
            (context / name).write_text(
                package.joinpath(name).read_text(encoding="utf-8"), encoding="utf-8"
            )
        run_command(
            ["docker", "build", "--tag", tag, str(context)],
            run_dir=run_dir,
            task="env.build.single-cell-python",
            backend="docker",
            inputs={
                "base": "python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b",
                "scanpy": "1.12.2",
            },
            outputs={"image": tag},
            artifacts=[
                {
                    "kind": "container-image",
                    "label": "single-cell-python",
                    "path": tag,
                }
            ],
        )
