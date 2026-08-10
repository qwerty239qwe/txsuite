from __future__ import annotations

import csv
import json
import math
import re
import tempfile
from importlib import resources
from pathlib import Path
from typing import Any

from txsuite.bulk import differential_expression_command
from txsuite.runtime import TxSuiteError, run_command

REQUIRED_COLUMNS = ("sample", "fastq_1", "fastq_2")
ALIGNERS = ("simpleaf", "star", "cellranger")
COMPARISON_COLUMNS = {
    "comparison",
    "group_column",
    "group_value",
    "design",
    "reference",
    "test",
    "method",
    "covariates",
    "padj",
    "lfc",
    "top_genes",
}
COMPARISON_INDEX_COLUMNS = (
    "comparison",
    "status",
    "method",
    "group_column",
    "group_value",
    "design",
    "reference",
    "test",
    "result",
    "significant",
    "error",
)


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
    check_inputs: bool = True,
    container_target: str | None = None,
    metadata: Path | None = None,
    barcode_column: str = "barcode",
    batch_column: str | None = None,
    integration: str = "none",
    counts_layer: str = "counts",
    target_sum: float = 10_000,
    n_hvg: int = 2_000,
    hvg_flavor: str = "seurat",
    n_pcs: int = 50,
    n_neighbors: int = 15,
    umap_min_dist: float = 0.5,
    marker_method: str = "wilcoxon",
    stop_after: str = "all",
    skip_umap: bool = False,
    skip_markers: bool = False,
    doublets: str = "off",
    doublet_batch_column: str | None = None,
    expected_doublet_rate: float = 0.05,
    doublet_threshold: float | None = None,
    top_markers: int = 100,
) -> list[str]:
    if check_inputs and not input_path.exists():
        raise TxSuiteError(f"Single-cell input does not exist: {input_path}")
    if not image.strip():
        raise TxSuiteError("Single-cell image cannot be empty")
    if min_genes < 0 or min_cells < 0 or not 0 <= max_mito_pct <= 100:
        raise TxSuiteError(
            "QC thresholds must be non-negative; max mito percent <= 100"
        )
    if resolution <= 0:
        raise TxSuiteError("Leiden resolution must be positive")
    if container_target is not None:
        if not re.fullmatch(r"/input/[A-Za-z0-9._-]+", container_target):
            raise TxSuiteError(
                "Single-cell container target must be a simple absolute path below /input"
            )
        target = container_target
    else:
        target = "/input/data" if input_path.is_dir() else f"/input/{input_path.name}"
    if metadata is not None and check_inputs and not metadata.is_file():
        raise TxSuiteError(f"Single-cell metadata does not exist: {metadata}")
    for label, value in (
        ("Barcode column", barcode_column),
        ("Batch column", batch_column),
        ("Doublet batch column", doublet_batch_column),
        ("Counts layer", counts_layer),
    ):
        if value is not None and not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", value):
            raise TxSuiteError(f"{label} must be a simple column name")
    if doublets not in {"off", "score", "filter"}:
        raise TxSuiteError("Doublet mode must be off, score, or filter")
    if integration not in {"none", "harmony"}:
        raise TxSuiteError("Integration method must be none or harmony")
    if integration == "harmony" and batch_column is None:
        raise TxSuiteError("Harmony integration requires a batch column")
    if integration == "harmony" and stop_after == "qc":
        raise TxSuiteError("Harmony integration requires PCA; stop-after cannot be qc")
    if not math.isfinite(target_sum) or target_sum <= 0:
        raise TxSuiteError("Normalization target sum must be positive")
    if n_hvg < 3 or n_pcs < 2 or n_neighbors < 2:
        raise TxSuiteError("HVGs must be >= 3; PCs and neighbors must be >= 2")
    if hvg_flavor not in {"seurat", "cell_ranger"}:
        raise TxSuiteError("HVG flavor must be seurat or cell_ranger")
    if not math.isfinite(umap_min_dist) or not 0 <= umap_min_dist <= 1:
        raise TxSuiteError("UMAP minimum distance must be in [0, 1]")
    if marker_method not in {"wilcoxon", "t-test"}:
        raise TxSuiteError("Marker method must be wilcoxon or t-test")
    if stop_after not in {"qc", "pca", "clusters", "all"}:
        raise TxSuiteError("Stop-after must be qc, pca, clusters, or all")
    if not math.isfinite(expected_doublet_rate) or not 0 < expected_doublet_rate < 1:
        raise TxSuiteError("Expected doublet rate must be in (0, 1)")
    if doublet_threshold is not None and (
        not math.isfinite(doublet_threshold) or doublet_threshold < 0
    ):
        raise TxSuiteError("Doublet threshold must be non-negative")
    if top_markers < 1:
        raise TxSuiteError("Top markers must be positive")
    command = [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={input_path.resolve()},target={target},readonly",
    ]
    if metadata is not None:
        command.extend(
            [
                "--mount",
                f"type=bind,source={metadata.resolve()},target=/input/metadata.tsv,readonly",
            ]
        )
    command.extend(
        [
            "--mount",
            f"type=bind,source={outdir.resolve()},target=/output",
            image,
            "python",
            "/opt/txsuite/single_cell.py",
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
            "--barcode-column",
            barcode_column,
            "--integration",
            integration,
            "--counts-layer",
            counts_layer,
            "--target-sum",
            str(target_sum),
            "--n-hvg",
            str(n_hvg),
            "--hvg-flavor",
            hvg_flavor,
            "--n-pcs",
            str(n_pcs),
            "--n-neighbors",
            str(n_neighbors),
            "--umap-min-dist",
            str(umap_min_dist),
            "--marker-method",
            marker_method,
            "--stop-after",
            stop_after,
            "--doublets",
            doublets,
            "--expected-doublet-rate",
            str(expected_doublet_rate),
            "--top-markers",
            str(top_markers),
        ]
    )
    if metadata is not None:
        command.extend(["--metadata", "/input/metadata.tsv"])
    if batch_column is not None:
        command.extend(["--batch-column", batch_column])
    if doublet_batch_column is not None:
        command.extend(["--doublet-batch-column", doublet_batch_column])
    if doublet_threshold is not None:
        command.extend(["--doublet-threshold", str(doublet_threshold)])
    if skip_umap:
        command.append("--skip-umap")
    if skip_markers:
        command.append("--skip-markers")
    return command


def pseudobulk_command(
    *,
    image: str,
    h5ad: Path,
    outdir: Path,
    sample_column: str,
    design: str,
    check_inputs: bool = True,
    group_column: str | None = None,
    group_value: str | None = None,
    covariates: tuple[str, ...] = (),
    reference: str | None = None,
    test: str | None = None,
) -> list[str]:
    if check_inputs and not h5ad.is_file():
        raise TxSuiteError(f"H5AD file does not exist: {h5ad}")
    if not image.strip():
        raise TxSuiteError("Single-cell image cannot be empty")
    for label, value in (
        ("Sample column", sample_column),
        ("Design", design),
        ("Group column", group_column),
        *(("Covariate", covariate) for covariate in covariates),
    ):
        if value is None:
            continue
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", value):
            raise TxSuiteError(f"{label} must be a simple column name")
    columns = (
        sample_column,
        design,
        *covariates,
        *([group_column] if group_column else []),
    )
    if len(set(columns)) != len(columns):
        raise TxSuiteError("Sample, design, group, and covariate columns must differ")
    if bool(group_column) != bool(group_value):
        raise TxSuiteError("Group column and group value must be used together")
    for label, value in (
        ("Reference", reference),
        ("Test", test),
        ("Group value", group_value),
    ):
        if value is not None and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
            raise TxSuiteError(f"{label} must be a simple factor level")
    if bool(reference) != bool(test) or reference is not None and reference == test:
        raise TxSuiteError("Reference and test must be used together and differ")
    command = [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={h5ad.resolve()},target=/input/data.h5ad,readonly",
        "--mount",
        f"type=bind,source={outdir.resolve()},target=/output",
        image,
        "python",
        "/opt/txsuite/single_cell.py",
        "pseudobulk",
        "/input/data.h5ad",
        "/output",
        "--sample-column",
        sample_column,
        "--design",
        design,
    ]
    if group_column is not None:
        command.extend(["--group-column", group_column, "--group-value", group_value])
    for covariate in covariates:
        command.extend(["--covariate", covariate])
    if reference is not None:
        command.extend(["--reference", reference, "--test", test])
    return command


def load_pseudobulk_manifest(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise TxSuiteError(f"Comparison manifest does not exist: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = set(reader.fieldnames or ())
        missing = {"comparison", "design", "reference", "test"} - headers
        unknown = headers - COMPARISON_COLUMNS
        if missing:
            raise TxSuiteError(
                f"Comparison manifest is missing columns: {', '.join(sorted(missing))}"
            )
        if unknown:
            raise TxSuiteError(
                f"Comparison manifest has unknown columns: {', '.join(sorted(unknown))}"
            )
        comparisons = []
        names = set()
        for line, raw in enumerate(reader, start=2):
            if None in raw:
                raise TxSuiteError(
                    f"Comparison manifest line {line} has more values than columns"
                )
            row = {key: (value or "").strip() for key, value in raw.items()}
            required = ("comparison", "design", "reference", "test")
            if any(not row[key] for key in required):
                raise TxSuiteError(
                    f"Comparison manifest line {line} has empty required values"
                )
            name = row["comparison"]
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
                raise TxSuiteError(f"Comparison manifest line {line} has invalid name")
            normalized_name = name.casefold()
            if normalized_name in names:
                raise TxSuiteError(f"Comparison name is duplicated: {name}")
            names.add(normalized_name)
            if bool(row.get("group_column")) != bool(row.get("group_value")):
                raise TxSuiteError(
                    f"Comparison manifest line {line} requires group column and value together"
                )
            method = row.get("method") or "deseq2"
            if method not in {"deseq2", "edger", "limma"}:
                raise TxSuiteError(
                    f"Comparison manifest line {line} has invalid method"
                )
            try:
                padj = float(row.get("padj") or 0.05)
                lfc = float(row.get("lfc") or 1)
                top_genes = int(row.get("top_genes") or 50)
            except ValueError as exc:
                raise TxSuiteError(
                    f"Comparison manifest line {line} has invalid numeric values"
                ) from exc
            comparisons.append(
                {
                    "comparison": name,
                    "group_column": row.get("group_column") or None,
                    "group_value": row.get("group_value") or None,
                    "design": row["design"],
                    "reference": row["reference"],
                    "test": row["test"],
                    "method": method,
                    "covariates": tuple(
                        value.strip()
                        for value in row.get("covariates", "").split(",")
                        if value.strip()
                    ),
                    "padj": padj,
                    "lfc": lfc,
                    "top_genes": top_genes,
                }
            )
    if not comparisons:
        raise TxSuiteError("Comparison manifest has no data rows")
    return comparisons


def _write_comparison_index(path: Path, records: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COMPARISON_INDEX_COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(records)


def _combine_comparison_results(records: list[dict[str, str]], output: Path) -> None:
    prefixes = [
        "comparison",
        "group_column",
        "group_value",
        "design",
        "reference",
        "test",
        "method",
    ]
    result_fields = []
    successful = [record for record in records if record["status"] != "failed"]
    for record in successful:
        with Path(record["result"]).open(encoding="utf-8", newline="") as handle:
            for field in csv.DictReader(handle, delimiter="\t").fieldnames or ():
                if field not in result_fields and field not in prefixes:
                    result_fields.append(field)
    with output.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(
            target, fieldnames=[*prefixes, *result_fields], delimiter="\t"
        )
        writer.writeheader()
        for record in successful:
            with Path(record["result"]).open(encoding="utf-8", newline="") as source:
                for row in csv.DictReader(source, delimiter="\t"):
                    writer.writerow({**{key: record[key] for key in prefixes}, **row})


def run_pseudobulk_manifest(
    *,
    manifest: Path,
    h5ad: Path,
    outdir: Path,
    sample_column: str,
    single_cell_image: str,
    bulk_image: str,
    resume: bool = False,
    dry_run: bool = False,
) -> list[list[str]]:
    comparisons = load_pseudobulk_manifest(manifest)
    commands = []
    records: list[dict[str, str]] = []
    if not dry_run:
        outdir.mkdir(parents=True, exist_ok=True)
        _write_comparison_index(outdir / "comparison-index.tsv", records)

    # ponytail: rows run sequentially; cache/group them when H5AD load time dominates.
    for comparison in comparisons:
        name = comparison["comparison"]
        comparison_out = outdir / name
        counts = comparison_out / "pseudobulk-counts.tsv"
        metadata = comparison_out / "pseudobulk-metadata.tsv"
        result = comparison_out / f"{comparison['method']}-results.tsv"
        significant = comparison_out / "significant-genes.tsv"
        run_root = comparison_out / ".txsuite"
        aggregate = pseudobulk_command(
            image=single_cell_image,
            h5ad=h5ad,
            outdir=comparison_out,
            sample_column=sample_column,
            design=comparison["design"],
            group_column=comparison["group_column"],
            group_value=comparison["group_value"],
            covariates=comparison["covariates"],
            reference=comparison["reference"],
            test=comparison["test"],
        )
        differential = differential_expression_command(
            method=comparison["method"],
            image=bulk_image,
            counts=counts,
            metadata=metadata,
            outdir=comparison_out,
            design=comparison["design"],
            reference=comparison["reference"],
            test=comparison["test"],
            covariates=comparison["covariates"],
            padj=comparison["padj"],
            lfc=comparison["lfc"],
            top_genes=comparison["top_genes"],
            check_inputs=False,
        )
        if dry_run:
            commands.extend([aggregate, differential])
            continue

        complete = False
        if resume and result.is_file() and significant.is_file():
            try:
                complete = (
                    json.loads((run_root / "de" / "run.json").read_text())["status"]
                    == "success"
                )
            except (OSError, KeyError, json.JSONDecodeError):
                pass
        status = "skipped" if complete else "success"
        error = ""
        try:
            if not complete:
                comparison_out.mkdir(parents=True, exist_ok=True)
                run_command(
                    aggregate,
                    run_dir=run_root / "pseudobulk",
                    task=f"single-cell.pseudobulk.{name}",
                    backend="Scanpy",
                    inputs={"manifest": str(manifest.resolve()), **comparison},
                    outputs={
                        "counts": str(counts.resolve()),
                        "metadata": str(metadata.resolve()),
                    },
                    artifacts=[
                        {
                            "kind": "table",
                            "label": "counts",
                            "path": str(counts.resolve()),
                        },
                        {
                            "kind": "table",
                            "label": "metadata",
                            "path": str(metadata.resolve()),
                        },
                    ],
                )
                if not counts.is_file() or not metadata.is_file():
                    raise TxSuiteError(
                        f"Comparison {name} did not produce pseudobulk tables"
                    )
                run_command(
                    differential,
                    run_dir=run_root / "de",
                    task=f"single-cell.pseudobulk-de.{name}",
                    backend={
                        "deseq2": "DESeq2",
                        "edger": "edgeR",
                        "limma": "limma-voom",
                    }[comparison["method"]],
                    inputs={
                        "manifest": str(manifest.resolve()),
                        "counts": str(counts.resolve()),
                        "metadata": str(metadata.resolve()),
                        **comparison,
                    },
                    outputs={"outdir": str(comparison_out.resolve())},
                    artifacts=[
                        {
                            "kind": "table",
                            "label": "DE results",
                            "path": str(result.resolve()),
                        },
                        {
                            "kind": "table",
                            "label": "significant genes",
                            "path": str(significant.resolve()),
                        },
                    ],
                )
                if not result.is_file() or not significant.is_file():
                    raise TxSuiteError(f"Comparison {name} did not produce DE tables")
        except TxSuiteError as exc:
            status = "failed"
            error = " ".join(str(exc).splitlines())
        records.append(
            {
                "comparison": name,
                "status": status,
                "method": comparison["method"],
                "group_column": comparison["group_column"] or "",
                "group_value": comparison["group_value"] or "",
                "design": comparison["design"],
                "reference": comparison["reference"],
                "test": comparison["test"],
                "result": str(result.resolve()) if status != "failed" else "",
                "significant": str(significant.resolve()) if status != "failed" else "",
                "error": error,
            }
        )
        _write_comparison_index(outdir / "comparison-index.tsv", records)

    if dry_run:
        return commands
    _combine_comparison_results(records, outdir / "combined-results.tsv")
    failed = sum(record["status"] == "failed" for record in records)
    if failed:
        raise TxSuiteError(
            f"{failed} comparisons failed; see {outdir / 'comparison-index.tsv'}"
        )
    return commands


def pseudobulk_workflow_command(
    config: dict[str, Any],
    *,
    h5ad: Path,
    outdir: Path,
    sample_column: str,
    design: str,
    reference: str,
    test: str,
    single_cell_image: str | None = None,
    bulk_image: str | None = None,
    padj: float = 0.05,
    lfc: float = 1.0,
    top_genes: int = 50,
    group_column: str | None = None,
    group_value: str | None = None,
    covariates: tuple[str, ...] = (),
    nextflow_config: Path | None = None,
    resume: bool = False,
    check_inputs: bool = True,
) -> list[str]:
    if check_inputs and not h5ad.is_file():
        raise TxSuiteError(f"H5AD file does not exist: {h5ad}")
    for label, value in (
        ("Sample column", sample_column),
        ("Design", design),
        ("Group column", group_column),
        *(("Covariate", covariate) for covariate in covariates),
    ):
        if value is None:
            continue
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", value):
            raise TxSuiteError(f"{label} must be a simple column name")
    columns = (
        sample_column,
        design,
        *covariates,
        *([group_column] if group_column else []),
    )
    if len(set(columns)) != len(columns):
        raise TxSuiteError("Sample, design, group, and covariate columns must differ")
    if bool(group_column) != bool(group_value):
        raise TxSuiteError("Group column and group value must be used together")
    for label, value in (
        ("Reference", reference),
        ("Test", test),
        ("Group value", group_value),
    ):
        if value is None:
            continue
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
            raise TxSuiteError(f"{label} must be a simple factor level")
    if reference == test:
        raise TxSuiteError("Reference and test levels must differ")
    if not 0 < padj <= 1 or lfc < 0 or top_genes < 1:
        raise TxSuiteError(
            "padj must be in (0, 1], lfc non-negative, and top genes positive"
        )
    command = _pseudobulk_nextflow_command(
        config,
        h5ad=h5ad,
        outdir=outdir,
        single_cell_image=single_cell_image,
        bulk_image=bulk_image,
        nextflow_config=nextflow_config,
        resume=resume,
        check_inputs=check_inputs,
    )
    command.extend(
        [
            "--sample_column",
            sample_column,
            "--design",
            design,
            "--reference",
            reference,
            "--test",
            test,
            "--padj",
            str(padj),
            "--lfc",
            str(lfc),
            "--top_genes",
            str(top_genes),
            "--covariates",
            ",".join(covariates),
        ]
    )
    if group_column is not None:
        command.extend(["--group_column", group_column, "--group_value", group_value])
    return command


def _pseudobulk_nextflow_command(
    config: dict[str, Any],
    *,
    h5ad: Path,
    outdir: Path,
    single_cell_image: str | None,
    bulk_image: str | None,
    nextflow_config: Path | None,
    resume: bool,
    check_inputs: bool = True,
) -> list[str]:
    if check_inputs and not h5ad.is_file():
        raise TxSuiteError(f"H5AD file does not exist: {h5ad}")
    images = config["images"]
    single_cell_image = single_cell_image or images["single_cell_python"]
    bulk_image = bulk_image or images["bulk_r"]
    if not single_cell_image.strip() or not bulk_image.strip():
        raise TxSuiteError("Workflow images cannot be empty")
    if nextflow_config is not None and not nextflow_config.is_file():
        raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")
    workflow = resources.files("txsuite.resources.nextflow").joinpath("main.nf")
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
            "--input",
            str(h5ad.resolve()),
            "--outdir",
            str(outdir.resolve()),
            "--single_cell_image",
            single_cell_image,
            "--bulk_image",
            bulk_image,
        ]
    )
    return command


def pseudobulk_manifest_workflow_command(
    config: dict[str, Any],
    *,
    manifest: Path,
    h5ad: Path,
    outdir: Path,
    sample_column: str,
    single_cell_image: str | None = None,
    bulk_image: str | None = None,
    nextflow_config: Path | None = None,
    resume: bool = False,
) -> list[str]:
    load_pseudobulk_manifest(manifest)
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", sample_column):
        raise TxSuiteError("Sample column must be a simple column name")
    command = _pseudobulk_nextflow_command(
        config,
        h5ad=h5ad,
        outdir=outdir,
        single_cell_image=single_cell_image,
        bulk_image=bulk_image,
        nextflow_config=nextflow_config,
        resume=resume,
    )
    command.extend(
        [
            "--manifest",
            str(manifest.resolve()),
            "--sample_column",
            sample_column,
        ]
    )
    return command


def cellranger_workflow_command(
    config: dict[str, Any],
    *,
    genome_name: str,
    fastqs: Path,
    sample: str,
    outdir: Path,
    fasta: Path | None = None,
    gtf: Path | None = None,
    reference: Path | None = None,
    cellranger_image: str | None = None,
    threads: int = 4,
    memory_gb: int = 8,
    create_bam: bool = True,
    nextflow_config: Path | None = None,
    resume: bool = False,
    check_inputs: bool = True,
) -> list[str]:
    """Build the native mkref+count Cell Ranger Nextflow DAG.

    Cell Ranger itself is never bundled: ``cellranger_image`` must name an image
    the caller already built with :func:`build_cellranger_image` from their own
    licensed download, and bare ``local``/``docker`` execution still requires a
    user-installed ``cellranger`` on ``PATH`` or inside that image.
    """

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", genome_name):
        raise TxSuiteError(
            "Cell Ranger genome name may contain only letters, numbers, ., _ and -"
        )
    if not re.fullmatch(r"[A-Za-z0-9_-]+", sample):
        raise TxSuiteError("Cell Ranger sample may contain only letters, numbers, _ and -")
    if bool(reference) == bool(fasta or gtf):
        raise TxSuiteError(
            "Pass either --reference, or both --fasta and --gtf to build one"
        )
    if bool(fasta) != bool(gtf):
        raise TxSuiteError("Building a reference requires both --fasta and --gtf")
    if threads < 1 or memory_gb < 1:
        raise TxSuiteError("Cell Ranger threads and memory must be positive")
    if check_inputs:
        if not fastqs.is_dir():
            raise TxSuiteError(f"FASTQ directory does not exist: {fastqs}")
        if reference is not None and not reference.is_dir():
            raise TxSuiteError(f"Cell Ranger reference does not exist: {reference}")
        if fasta is not None and not fasta.is_file():
            raise TxSuiteError(f"Genome FASTA does not exist: {fasta}")
        if gtf is not None and not gtf.is_file():
            raise TxSuiteError(f"Annotation GTF does not exist: {gtf}")
    if nextflow_config is not None and not nextflow_config.is_file():
        raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")

    workflow = resources.files("txsuite.resources.nextflow").joinpath(
        "single_cell_cellranger.nf"
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
            "--genome_name",
            genome_name,
            "--fastqs",
            str(fastqs.resolve()),
            "--sample",
            sample,
            "--outdir",
            str(outdir.resolve()),
            "--cellranger_threads",
            str(threads),
            "--cellranger_memory_gb",
            str(memory_gb),
            "--cellranger_create_bam",
            "true" if create_bam else "false",
        ]
    )
    if reference is not None:
        command.extend(["--cellranger_reference", str(reference.resolve())])
    else:
        command.extend(
            [
                "--fasta",
                str(fasta.resolve()),
                "--gtf",
                str(gtf.resolve()),
            ]
        )
    if cellranger_image:
        command.extend(["--cellranger_image", cellranger_image])
    return command


def build_cellranger_image(tag: str, *, source_tarball: Path, run_dir: Path) -> None:
    """Build a local Cell Ranger image from a tarball the caller already licensed.

    TxSuite packages only the install recipe. ``source_tarball`` must already be
    on disk, downloaded by the caller from their own 10x Genomics account; it is
    copied into the build context and never fetched or redistributed by TxSuite.
    """

    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    if not source_tarball.is_file():
        raise TxSuiteError(f"Cell Ranger tarball does not exist: {source_tarball}")
    package = resources.files("txsuite.resources.cellranger")
    with tempfile.TemporaryDirectory(prefix="txsuite-cellranger-") as directory:
        context = Path(directory)
        (context / "Dockerfile").write_text(
            package.joinpath("Dockerfile").read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        tarball_name = source_tarball.name
        with (context / tarball_name).open("wb") as destination:
            with source_tarball.open("rb") as source:
                destination.write(source.read())
        run_command(
            [
                "docker",
                "build",
                "--tag",
                tag,
                "--build-arg",
                f"CELLRANGER_TARBALL={tarball_name}",
                str(context),
            ],
            run_dir=run_dir,
            task="env.build.cellranger",
            backend="docker",
            inputs={"source_tarball": str(source_tarball.resolve())},
            outputs={"image": tag},
            artifacts=[
                {
                    "kind": "container-image",
                    "label": "cellranger",
                    "path": tag,
                }
            ],
        )


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
