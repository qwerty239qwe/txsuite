from __future__ import annotations

import re
import tempfile
from importlib import resources
from pathlib import Path

from txsuite.runtime import TxSuiteError, run_command


def spaceranger_command(
    *,
    run_id: str,
    transcriptome: Path,
    fastqs: Path,
    image: Path,
    sample: str | None,
    slide: str | None,
    area: str | None,
    unknown_slide: bool,
    create_bam: bool,
    cores: int | None,
    memory: int | None,
) -> list[str]:
    if not re.fullmatch(r"[A-Za-z0-9_-]+", run_id):
        raise TxSuiteError("Space Ranger id may contain only letters, numbers, _ and -")
    for label, path, kind in (
        ("Transcriptome", transcriptome, "directory"),
        ("FASTQ", fastqs, "directory"),
        ("Image", image, "file"),
    ):
        exists = path.is_dir() if kind == "directory" else path.is_file()
        if not exists:
            raise TxSuiteError(f"{label} {kind} does not exist: {path}")
    if unknown_slide == bool(slide or area) or bool(slide) != bool(area):
        raise TxSuiteError("Use --slide with --area, or use --unknown-slide")
    if cores is not None and cores < 1 or memory is not None and memory < 1:
        raise TxSuiteError("Space Ranger cores and memory must be positive")

    command = [
        "spaceranger",
        "count",
        f"--id={run_id}",
        f"--transcriptome={transcriptome.resolve()}",
        f"--fastqs={fastqs.resolve()}",
        f"--image={image.resolve()}",
        f"--create-bam={'true' if create_bam else 'false'}",
    ]
    if sample:
        command.append(f"--sample={sample}")
    if unknown_slide:
        command.append("--unknown-slide")
    else:
        command.extend((f"--slide={slide}", f"--area={area}"))
    if cores is not None:
        command.append(f"--localcores={cores}")
    if memory is not None:
        command.append(f"--localmem={memory}")
    return command


def spacemake_command(project_root: Path, cores: int) -> list[str]:
    if not project_root.is_dir():
        raise TxSuiteError(f"Spacemake project does not exist: {project_root}")
    if cores < 1:
        raise TxSuiteError("Spacemake cores must be positive")
    return ["spacemake", "run", "--cores", str(cores), "--keep-going"]


def analysis_command(
    *,
    image: str,
    input_path: Path,
    outdir: Path,
    dataset_id: str,
    min_counts: int,
    min_spots: int,
) -> list[str]:
    if not input_path.is_dir():
        raise TxSuiteError(f"Spatial input directory does not exist: {input_path}")
    if not image.strip():
        raise TxSuiteError("Spatial image cannot be empty")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", dataset_id):
        raise TxSuiteError("Dataset id may contain only letters, numbers, _ and -")
    if min_counts < 0 or min_spots < 0:
        raise TxSuiteError("Spatial QC thresholds must be non-negative")
    target = (
        "/input/data.zarr" if input_path.name.endswith(".zarr") else "/input/visium"
    )
    return [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={input_path.resolve()},target={target},readonly",
        "--mount",
        f"type=bind,source={outdir.resolve()},target=/output",
        image,
        target,
        "/output",
        "--dataset-id",
        dataset_id,
        "--min-counts",
        str(min_counts),
        "--min-spots",
        str(min_spots),
    ]


def build_spatial_image(tag: str, *, run_dir: Path) -> None:
    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    package = resources.files("txsuite.resources.spatial_python")
    with tempfile.TemporaryDirectory(prefix="txsuite-spatial-") as directory:
        context = Path(directory)
        for name in ("Dockerfile", "spatial.py"):
            (context / name).write_text(
                package.joinpath(name).read_text(encoding="utf-8"), encoding="utf-8"
            )
        run_command(
            ["docker", "build", "--tag", tag, str(context)],
            run_dir=run_dir,
            task="env.build.spatial-python",
            backend="docker",
            inputs={
                "base": "python:3.12-slim-bookworm@sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b",
                "spatialdata": "0.8.0",
                "spatialdata-io": "0.7.1",
                "squidpy": "1.8.3",
            },
            outputs={"image": tag},
            artifacts=[
                {
                    "kind": "container-image",
                    "label": "spatial-python",
                    "path": tag,
                }
            ],
        )
