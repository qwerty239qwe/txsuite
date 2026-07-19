from __future__ import annotations

import csv
import math
import re
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError, run_command


REQUIRED_COLUMNS = {"sample", "fastq_1", "fastq_2", "strandedness"}
STRANDEDNESS = {"auto", "forward", "reverse", "unstranded"}
DE_METHODS = ("deseq2", "edger", "limma")


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
    covariates: tuple[str, ...] = (),
    padj: float = 0.05,
    lfc: float = 1.0,
    top_genes: int = 50,
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
    if any(
        not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", covariate)
        for covariate in covariates
    ):
        raise TxSuiteError("Covariates must be simple metadata column names")
    if design in covariates or len(set(covariates)) != len(covariates):
        raise TxSuiteError("Covariates must be unique and different from design")
    if not math.isfinite(padj) or not 0 < padj <= 1:
        raise TxSuiteError("Adjusted p-value threshold must be in (0, 1]")
    if not math.isfinite(lfc) or lfc < 0:
        raise TxSuiteError("Absolute log2 fold-change threshold must be non-negative")
    if top_genes < 1:
        raise TxSuiteError("Top genes must be positive")
    if not image.strip():
        raise TxSuiteError("Bulk R image cannot be empty")
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
        "Rscript",
        "/opt/txsuite/deseq2.R",
        "/input/counts.tsv",
        "/input/metadata.tsv",
        design,
        reference,
        test,
        "/output",
        str(padj),
        str(lfc),
        str(top_genes),
        ",".join(covariates),
    ]


def differential_expression_command(
    *,
    method: str,
    image: str,
    counts: Path,
    metadata: Path,
    design: str,
    reference: str,
    test: str,
    outdir: Path,
    covariates: tuple[str, ...] = (),
    padj: float = 0.05,
    lfc: float = 1.0,
    top_genes: int = 50,
    check_inputs: bool = True,
) -> list[str]:
    if method not in DE_METHODS:
        raise TxSuiteError(f"DE method must be one of: {', '.join(DE_METHODS)}")
    command = deseq2_command(
        image=image,
        counts=counts,
        metadata=metadata,
        design=design,
        reference=reference,
        test=test,
        outdir=outdir,
        covariates=covariates,
        padj=padj,
        lfc=lfc,
        top_genes=top_genes,
        check_inputs=check_inputs,
    )
    if method != "deseq2":
        script = command.index("/opt/txsuite/deseq2.R")
        command[script] = "/opt/txsuite/alternative_de.R"
        command.insert(script + 1, method)
    return command


def enrichment_command(
    *,
    image: str,
    de_results: Path,
    genesets: Path,
    mode: str,
    outdir: Path,
    padj: float = 0.05,
    lfc: float = 1.0,
    min_size: int = 10,
    max_size: int = 500,
    adjust: str = "BH",
) -> list[str]:
    for label, path in (("DE results", de_results), ("GMT gene sets", genesets)):
        if not path.is_file():
            raise TxSuiteError(f"{label} file does not exist: {path}")
    if mode not in {"ora", "gsea"}:
        raise TxSuiteError("Enrichment mode must be 'ora' or 'gsea'")
    if not image.strip():
        raise TxSuiteError("Bulk R image cannot be empty")
    if not math.isfinite(padj) or not 0 < padj <= 1:
        raise TxSuiteError("Adjusted p-value threshold must be in (0, 1]")
    if not math.isfinite(lfc) or lfc < 0:
        raise TxSuiteError("Absolute log2 fold-change threshold must be non-negative")
    if min_size < 1 or max_size < min_size:
        raise TxSuiteError("Gene-set sizes must satisfy 1 <= min-size <= max-size")
    if adjust not in {
        "holm",
        "hochberg",
        "hommel",
        "bonferroni",
        "BH",
        "BY",
        "fdr",
        "none",
    }:
        raise TxSuiteError(f"Unknown p-value adjustment method: {adjust}")
    mounts = (
        f"type=bind,source={de_results.resolve()},target=/input/de.tsv,readonly",
        f"type=bind,source={genesets.resolve()},target=/input/genesets.gmt,readonly",
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
        "Rscript",
        "/opt/txsuite/enrichment.R",
        mode,
        "/input/de.tsv",
        "/input/genesets.gmt",
        "/output",
        str(padj),
        str(lfc),
        str(min_size),
        str(max_size),
        adjust,
    ]


def build_bulk_r_image(tag: str, *, run_dir: Path) -> None:
    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    package = resources.files("txsuite.resources.bulk_r")
    with tempfile.TemporaryDirectory(prefix="txsuite-bulk-r-") as directory:
        context = Path(directory)
        for name in ("Dockerfile", "deseq2.R", "alternative_de.R", "enrichment.R"):
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
                ),
                "DESeq2": "1.52.0",
                "edgeR": "4.10.1",
                "limma": "3.68.4",
                "clusterProfiler": "4.20.0",
            },
            outputs={"image": tag},
            artifacts=[{"kind": "container-image", "label": "bulk-r", "path": tag}],
        )
