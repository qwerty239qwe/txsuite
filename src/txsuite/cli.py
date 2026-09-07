from __future__ import annotations

import argparse
import json
import shlex
import shutil
from dataclasses import replace
from pathlib import Path

from txsuite.bulk import (
    build_bulk_r_image,
    deseq2_command,
    differential_expression_command,
    enrichment_command,
)
from txsuite.bulk import (
    validate_samplesheet as validate_bulk_samplesheet,
)
from txsuite.bulk import (
    workflow_command as bulk_workflow_command,
)
from txsuite.catalog import select_tools
from txsuite.config import DEFAULT_TOML, ConfigError, load_config
from txsuite.hardening import cache_reference, image_is_locked
from txsuite.project import (
    ProjectExecutor,
    ProvenanceError,
    RunBundle,
    StageExecutionError,
    WorkflowConfigError,
    hash_payload,
    load_project_config,
    plan_workflow,
    render_plan,
    scaffold_project_preset,
    select_stages,
)
from txsuite.runtime import TxSuiteError, format_command, run_command
from txsuite.single_cell import (
    analysis_command,
    build_cellranger_image,
    build_single_cell_image,
    cellranger_workflow_command,
    pseudobulk_command,
    pseudobulk_manifest_workflow_command,
    pseudobulk_workflow_command,
    run_pseudobulk_manifest,
)
from txsuite.single_cell import (
    validate_samplesheet as validate_single_cell_samplesheet,
)
from txsuite.single_cell import (
    workflow_command as single_cell_workflow_command,
)
from txsuite.spatial import (
    analysis_command as spatial_analysis_command,
)
from txsuite.spatial import (
    build_spatial_image,
    spacemake_command,
    spaceranger_command,
)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="txsuite", description="Transcriptomics toolbox"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    tools = commands.add_parser("tools", help="list selected analysis tools")
    tools.add_argument("--modality", choices=("bulk", "single-cell", "spatial"))
    tools.add_argument(
        "--stage",
        choices=(
            "workflow",
            "qc",
            "trim",
            "alignment",
            "quantification",
            "differential-expression",
        ),
    )

    config = commands.add_parser("config", help="manage layered TOML configuration")
    config_commands = config.add_subparsers(dest="config_command", required=True)
    show = config_commands.add_parser("show", help="show the merged configuration")
    show.add_argument("--file", type=Path, default=Path("txsuite.toml"))
    init = config_commands.add_parser("init", help="write a project configuration")
    init.add_argument("--file", type=Path, default=Path("txsuite.toml"))
    init.add_argument("--force", action="store_true")

    project = commands.add_parser(
        "project", help="scaffold, inspect, and run project workflows"
    )
    project_commands = project.add_subparsers(
        dest="project_command", required=True
    )
    project_init = project_commands.add_parser(
        "init", help="scaffold a packaged project preset"
    )
    project_init.add_argument("--preset", required=True)
    project_init.add_argument("target", type=Path)

    project_validate = project_commands.add_parser(
        "validate", help="validate and fully plan a project workflow"
    )
    project_validate.add_argument("workflow", type=Path)
    project_validate.add_argument(
        "--config", type=Path, default=Path("txsuite.toml")
    )

    project_plan = project_commands.add_parser(
        "plan", help="render a side-effect-free project command plan"
    )
    project_plan.add_argument("workflow", type=Path)
    project_plan.add_argument(
        "--config", type=Path, default=Path("txsuite.toml")
    )
    project_plan.add_argument("--json", action="store_true")

    project_run = project_commands.add_parser(
        "run", help="plan and execute a project workflow"
    )
    project_run.add_argument("workflow", type=Path)
    project_run.add_argument(
        "--config", type=Path, default=Path("txsuite.toml")
    )
    project_run.add_argument("--dry-run", action="store_true")
    project_run.add_argument("--resume", action="store_true")
    project_run.add_argument("--from", dest="from_stage")
    project_run.add_argument("--to", dest="to_stage")
    project_run.add_argument("--stages")
    project_run.add_argument("--run-id")

    project_status = project_commands.add_parser(
        "status", help="show the immutable manifest view of a project run"
    )
    project_status.add_argument("run_dir", type=Path)

    workflow = commands.add_parser(
        "workflow", help="run complete transcriptomics workflows"
    )
    workflow_commands = workflow.add_subparsers(dest="workflow_command", required=True)
    bulk_workflow = workflow_commands.add_parser(
        "bulk", help="run pinned nf-core/rnaseq"
    )
    bulk_workflow.add_argument("--input", type=Path, required=True)
    bulk_workflow.add_argument("--outdir", type=Path, required=True)
    bulk_workflow.add_argument("--params-file", type=Path)
    bulk_workflow.add_argument("--nextflow-config", type=Path)
    bulk_workflow.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    bulk_workflow.add_argument("--run-dir", type=Path)
    bulk_workflow.add_argument("--resume", action="store_true")
    bulk_workflow.add_argument("--dry-run", action="store_true")

    single_cell_workflow = workflow_commands.add_parser(
        "single-cell", help="run pinned nf-core/scrnaseq"
    )
    single_cell_workflow.add_argument("--input", type=Path, required=True)
    single_cell_workflow.add_argument("--outdir", type=Path, required=True)
    single_cell_workflow.add_argument(
        "--aligner", choices=("simpleaf", "star", "cellranger"), default="simpleaf"
    )
    single_cell_workflow.add_argument("--protocol")
    single_cell_workflow.add_argument("--params-file", type=Path)
    single_cell_workflow.add_argument("--nextflow-config", type=Path)
    single_cell_workflow.add_argument(
        "--config", type=Path, default=Path("txsuite.toml")
    )
    single_cell_workflow.add_argument("--run-dir", type=Path)
    single_cell_workflow.add_argument("--resume", action="store_true")
    single_cell_workflow.add_argument("--dry-run", action="store_true")

    pseudobulk_workflow = workflow_commands.add_parser(
        "pseudobulk-de", help="run pseudobulk and DESeq2 as a native Nextflow DAG"
    )
    pseudobulk_workflow.add_argument("--input", type=Path, required=True)
    pseudobulk_workflow.add_argument("--sample-column", required=True)
    pseudobulk_workflow.add_argument("--counts-layer", default="counts")
    pseudobulk_workflow.add_argument("--design", required=True)
    pseudobulk_workflow.add_argument("--group-column")
    pseudobulk_workflow.add_argument("--group-value")
    pseudobulk_workflow.add_argument("--covariate", action="append", default=[])
    pseudobulk_workflow.add_argument("--reference", required=True)
    pseudobulk_workflow.add_argument("--test", required=True)
    pseudobulk_workflow.add_argument("--padj", type=float, default=0.05)
    pseudobulk_workflow.add_argument("--lfc", type=float, default=1.0)
    pseudobulk_workflow.add_argument("--top-genes", type=int, default=50)
    pseudobulk_workflow.add_argument("--outdir", type=Path, required=True)
    pseudobulk_workflow.add_argument("--image")
    pseudobulk_workflow.add_argument("--bulk-image")
    pseudobulk_workflow.add_argument("--nextflow-config", type=Path)
    pseudobulk_workflow.add_argument(
        "--config", type=Path, default=Path("txsuite.toml")
    )
    pseudobulk_workflow.add_argument("--run-dir", type=Path)
    pseudobulk_workflow.add_argument("--resume", action="store_true")
    pseudobulk_workflow.add_argument("--dry-run", action="store_true")

    cellranger_workflow = workflow_commands.add_parser(
        "single-cell-cellranger",
        help="run user-installed Cell Ranger as a native mkref+count Nextflow DAG",
    )
    cellranger_workflow.add_argument("--genome-name", required=True)
    cellranger_workflow.add_argument("--fastqs", type=Path, required=True)
    cellranger_workflow.add_argument("--sample", required=True)
    cellranger_workflow.add_argument("--fasta", type=Path)
    cellranger_workflow.add_argument("--gtf", type=Path)
    cellranger_workflow.add_argument("--reference", type=Path)
    cellranger_workflow.add_argument("--cellranger-image")
    cellranger_workflow.add_argument("--threads", type=int, default=4)
    cellranger_workflow.add_argument("--memory-gb", type=int, default=8)
    cellranger_workflow.add_argument("--no-bam", action="store_true")
    cellranger_workflow.add_argument("--outdir", type=Path, required=True)
    cellranger_workflow.add_argument("--nextflow-config", type=Path)
    cellranger_workflow.add_argument(
        "--config", type=Path, default=Path("txsuite.toml")
    )
    cellranger_workflow.add_argument("--run-dir", type=Path)
    cellranger_workflow.add_argument("--resume", action="store_true")
    cellranger_workflow.add_argument("--dry-run", action="store_true")

    spatial_workflow = workflow_commands.add_parser(
        "spatial", help="run user-installed Space Ranger 4.1"
    )
    spatial_workflow.add_argument("--id", required=True)
    spatial_workflow.add_argument("--transcriptome", type=Path, required=True)
    spatial_workflow.add_argument("--fastqs", type=Path, required=True)
    spatial_workflow.add_argument("--image", type=Path, required=True)
    spatial_workflow.add_argument("--sample")
    spatial_workflow.add_argument("--slide")
    spatial_workflow.add_argument("--area")
    spatial_workflow.add_argument("--unknown-slide", action="store_true")
    spatial_workflow.add_argument("--no-bam", action="store_true")
    spatial_workflow.add_argument("--cores", type=int)
    spatial_workflow.add_argument("--memory", type=int)
    spatial_workflow.add_argument("--outdir", type=Path, required=True)
    spatial_workflow.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    spatial_workflow.add_argument("--run-dir", type=Path)
    spatial_workflow.add_argument("--dry-run", action="store_true")

    spatial_open = workflow_commands.add_parser(
        "spatial-open", help="run an existing Spacemake project (experimental)"
    )
    spatial_open.add_argument("--project-root", type=Path, required=True)
    spatial_open.add_argument("--cores", type=int, default=1)
    spatial_open.add_argument("--run-dir", type=Path)
    spatial_open.add_argument("--dry-run", action="store_true")

    bulk = commands.add_parser("bulk", help="bulk RNA-seq downstream analysis")
    bulk_commands = bulk.add_subparsers(dest="bulk_command", required=True)
    de = bulk_commands.add_parser("de", help="run differential expression")
    de.add_argument("--method", choices=("deseq2", "edger", "limma"), default="deseq2")
    de.add_argument("--counts", type=Path, required=True)
    de.add_argument("--metadata", type=Path, required=True)
    de.add_argument("--design")
    de.add_argument("--reference")
    de.add_argument("--test")
    de.add_argument("--covariate", action="append", default=[])
    de.add_argument("--formula")
    de.add_argument("--coefficient")
    de.add_argument("--padj", type=float, default=0.05)
    de.add_argument("--lfc", type=float, default=1.0)
    de.add_argument("--top-genes", type=int, default=50)
    de.add_argument("--outdir", type=Path, required=True)
    de.add_argument("--image")
    de.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    de.add_argument("--run-dir", type=Path)
    de.add_argument("--dry-run", action="store_true")

    enrich = bulk_commands.add_parser(
        "enrich", help="run GMT-based over-representation analysis or GSEA"
    )
    enrich.add_argument("--de", type=Path, required=True)
    enrich.add_argument("--genesets", type=Path, required=True)
    enrich.add_argument("--mode", choices=("ora", "gsea"), required=True)
    enrich.add_argument("--padj", type=float, default=0.05)
    enrich.add_argument("--lfc", type=float, default=1.0)
    enrich.add_argument("--min-size", type=int, default=10)
    enrich.add_argument("--max-size", type=int, default=500)
    enrich.add_argument(
        "--adjust",
        choices=("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
        default="BH",
    )
    enrich.add_argument("--outdir", type=Path, required=True)
    enrich.add_argument("--image")
    enrich.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    enrich.add_argument("--run-dir", type=Path)
    enrich.add_argument("--dry-run", action="store_true")

    single_cell = commands.add_parser(
        "single-cell", help="single-cell RNA-seq downstream analysis"
    )
    single_cell_commands = single_cell.add_subparsers(
        dest="single_cell_command", required=True
    )
    analyze = single_cell_commands.add_parser(
        "analyze", help="run Scanpy QC, normalization, clustering, and UMAP"
    )
    analyze.add_argument("--input", type=Path, required=True)
    analyze.add_argument("--outdir", type=Path, required=True)
    analyze.add_argument("--min-genes", type=int, default=200)
    analyze.add_argument("--min-cells", type=int, default=3)
    analyze.add_argument("--max-mito-pct", type=float, default=20)
    analyze.add_argument("--resolution", type=float, default=1)
    analyze.add_argument("--metadata", type=Path)
    analyze.add_argument("--barcode-column", default="barcode")
    analyze.add_argument("--batch-column")
    analyze.add_argument("--integration", choices=("none", "harmony"), default="none")
    analyze.add_argument("--counts-layer", default="counts")
    analyze.add_argument("--target-sum", type=float, default=10_000)
    analyze.add_argument("--n-hvg", type=int, default=2_000)
    analyze.add_argument(
        "--hvg-flavor", choices=("seurat", "cell_ranger"), default="seurat"
    )
    analyze.add_argument("--n-pcs", type=int, default=50)
    analyze.add_argument("--n-neighbors", type=int, default=15)
    analyze.add_argument("--umap-min-dist", type=float, default=0.5)
    analyze.add_argument(
        "--marker-method", choices=("wilcoxon", "t-test"), default="wilcoxon"
    )
    analyze.add_argument(
        "--stop-after", choices=("qc", "pca", "clusters", "all"), default="all"
    )
    analyze.add_argument("--skip-umap", action="store_true")
    analyze.add_argument("--skip-markers", action="store_true")
    analyze.add_argument(
        "--doublets", choices=("off", "score", "filter"), default="off"
    )
    analyze.add_argument("--doublet-batch-column")
    analyze.add_argument("--expected-doublet-rate", type=float, default=0.05)
    analyze.add_argument("--doublet-threshold", type=float)
    analyze.add_argument("--top-markers", type=int, default=100)
    analyze.add_argument("--image")
    analyze.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    analyze.add_argument("--run-dir", type=Path)
    analyze.add_argument("--dry-run", action="store_true")

    pseudobulk = single_cell_commands.add_parser(
        "pseudobulk-de", help="aggregate cells and run the bulk DESeq2 backend"
    )
    pseudobulk.add_argument("--input", type=Path, required=True)
    pseudobulk.add_argument("--sample-column", required=True)
    pseudobulk.add_argument("--counts-layer", default="counts")
    pseudobulk.add_argument("--design", required=True)
    pseudobulk.add_argument("--group-column")
    pseudobulk.add_argument("--group-value")
    pseudobulk.add_argument("--covariate", action="append", default=[])
    pseudobulk.add_argument("--reference", required=True)
    pseudobulk.add_argument("--test", required=True)
    pseudobulk.add_argument("--outdir", type=Path, required=True)
    pseudobulk.add_argument("--image")
    pseudobulk.add_argument("--bulk-image")
    pseudobulk.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    pseudobulk.add_argument("--run-dir", type=Path)
    pseudobulk.add_argument("--dry-run", action="store_true")

    pseudobulk_batch = single_cell_commands.add_parser(
        "pseudobulk-batch",
        help="run groups and contrasts from a comparison manifest",
    )
    pseudobulk_batch.add_argument("--input", type=Path, required=True)
    pseudobulk_batch.add_argument("--sample-column", required=True)
    pseudobulk_batch.add_argument("--counts-layer", default="counts")
    pseudobulk_batch.add_argument("--manifest", type=Path, required=True)
    pseudobulk_batch.add_argument("--outdir", type=Path, required=True)
    pseudobulk_batch.add_argument("--image")
    pseudobulk_batch.add_argument("--bulk-image")
    pseudobulk_batch.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    pseudobulk_batch.add_argument("--nextflow-config", type=Path)
    pseudobulk_batch.add_argument("--run-dir", type=Path)
    pseudobulk_batch.add_argument(
        "--direct", action="store_true", help="run each comparison directly with Docker"
    )
    pseudobulk_batch.add_argument("--resume", action="store_true")
    pseudobulk_batch.add_argument("--dry-run", action="store_true")

    spatial = commands.add_parser("spatial", help="spatial transcriptomics analysis")
    spatial_commands = spatial.add_subparsers(dest="spatial_command", required=True)
    spatial_analyze = spatial_commands.add_parser(
        "analyze", help="import Visium and build a Squidpy spatial graph"
    )
    spatial_analyze.add_argument("--input", type=Path, required=True)
    spatial_analyze.add_argument("--outdir", type=Path, required=True)
    spatial_analyze.add_argument("--dataset-id", default="sample")
    spatial_analyze.add_argument("--min-counts", type=int, default=500)
    spatial_analyze.add_argument("--min-spots", type=int, default=3)
    spatial_analyze.add_argument("--image")
    spatial_analyze.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    spatial_analyze.add_argument("--run-dir", type=Path)
    spatial_analyze.add_argument("--dry-run", action="store_true")

    reference = commands.add_parser("reference", help="manage verified references")
    reference_commands = reference.add_subparsers(
        dest="reference_command", required=True
    )
    reference_cache = reference_commands.add_parser(
        "cache", help="cache a file only when its SHA-256 matches"
    )
    reference_cache.add_argument("--source", required=True)
    reference_cache.add_argument("--sha256", required=True)
    reference_cache.add_argument("--name", required=True)
    reference_cache.add_argument(
        "--root", type=Path, default=Path(".txsuite/references")
    )

    env = commands.add_parser("env", help="inspect execution environments")
    env_commands = env.add_subparsers(dest="env_command", required=True)
    list_environments = env_commands.add_parser(
        "list", help="list configured TxSuite-owned images"
    )
    list_environments.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    doctor = env_commands.add_parser("doctor", help="check workflow prerequisites")
    doctor.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    verify_images = env_commands.add_parser(
        "verify-images", help="fail when configured images use mutable tags"
    )
    verify_images.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    build = env_commands.add_parser("build", help="build a TxSuite-owned image")
    build.add_argument(
        "environment",
        choices=("bulk-r", "single-cell-python", "spatial-python", "cellranger"),
    )
    build.add_argument("--tag")
    build.add_argument(
        "--source-tarball",
        type=Path,
        help="user-downloaded cellranger-*.tar.gz (required for 'cellranger')",
    )
    build.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    build.add_argument("--run-dir", type=Path)
    build.add_argument("--dry-run", action="store_true")
    return parser


def _load_project_plan(workflow_path: Path, config_path: Path):
    """Load both configuration layers and cross the full planning boundary."""

    workflow = load_project_config(workflow_path)
    global_config = load_config(config_path)
    return workflow, global_config, plan_workflow(workflow, global_config)


def _selected_project_stages(value: str | None) -> tuple[str, ...] | None:
    if value is None:
        return None
    selected = tuple(item.strip() for item in value.split(",") if item.strip())
    if not selected:
        raise TxSuiteError("--stages must contain at least one stage ID")
    if len(selected) != len(set(selected)):
        raise TxSuiteError("--stages must not contain duplicate stage IDs")
    return selected


def _effective_project_stages(
    plan,
    *,
    selected_stage_ids: tuple[str, ...] | None,
    from_stage: str | None,
    to_stage: str | None,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return requested stages plus their transitive prerequisites in plan order."""

    stages = plan.to_dict()["stages"]
    requested = select_stages(
        stages,
        selected_stage_ids=selected_stage_ids,
        from_stage=from_stage,
        to_stage=to_stage,
    )
    requested_ids = tuple(str(stage["id"]) for stage in requested)
    dependencies = {
        str(stage["id"]): tuple(str(item) for item in stage.get("depends_on", ()))
        for stage in stages
    }
    included = set(requested_ids)
    pending = list(requested_ids)
    while pending:
        stage_id = pending.pop()
        for dependency in dependencies[stage_id]:
            if dependency not in included:
                included.add(dependency)
                pending.append(dependency)

    effective = tuple(
        str(stage["id"]) for stage in stages if str(stage["id"]) in included
    )
    requested_set = set(requested_ids)
    auto_added = tuple(stage_id for stage_id in effective if stage_id not in requested_set)
    return effective, auto_added


def _render_project_run_plan(
    plan, effective_stage_ids: tuple[str, ...], auto_added: tuple[str, ...]
) -> str:
    included = set(effective_stage_ids)
    selected_plan = replace(
        plan,
        stages=tuple(stage for stage in plan.stages if stage.id in included),
    )
    lines = [f"Selected stages: {', '.join(effective_stage_ids)}"]
    if auto_added:
        lines.append(f"Auto-added prerequisites: {', '.join(auto_added)}")
    lines.extend(("", render_plan(selected_plan, format="human")))
    return "\n".join(lines)


def _project_runs_dir(workflow) -> Path:
    return workflow.project.output_root / workflow.project.id / "runs"


def _load_project_bundle(run_dir: Path) -> tuple[RunBundle, dict]:
    try:
        bundle = RunBundle.load(run_dir)
        manifest = bundle.manifest()
    except (OSError, ValueError, KeyError, ProvenanceError) as exc:
        raise TxSuiteError(f"Cannot read run bundle {run_dir}: {exc}") from exc
    if (
        not isinstance(manifest.get("run_id"), str)
        or not isinstance(manifest.get("status"), str)
        or not isinstance(manifest.get("hashes"), dict)
    ):
        raise TxSuiteError(
            f"Invalid run manifest {bundle.manifest_path}: "
            "run_id, status, and hashes are required"
        )
    return bundle, manifest


def _require_compatible_resume(
    bundle: RunBundle, manifest: dict, workflow, global_config: dict, plan
) -> None:
    stored = manifest["hashes"]
    current = {
        "workflow": hash_payload(workflow.to_dict()),
        "resolved_config": hash_payload(global_config),
        "command_plan": hash_payload(plan.to_dict()),
    }
    labels = {
        "workflow": "workflow",
        "resolved_config": "resolved configuration",
        "command_plan": "command plan",
    }
    changed = [
        labels[name]
        for name, digest in current.items()
        if stored.get(name) != digest
    ]
    if changed:
        raise TxSuiteError(
            f"Cannot resume run {bundle.run_id!r}: immutable "
            f"{', '.join(changed)} snapshot changed. "
            "Start a new run with a new --run-id."
        )


def _run_project_command(args: argparse.Namespace) -> int:
    if args.project_command == "init":
        print(scaffold_project_preset(args.preset, args.target))
        return 0

    if args.project_command == "status":
        _, manifest = _load_project_bundle(args.run_dir)
        print(json.dumps(manifest, indent=2, sort_keys=True))
        return 0

    workflow, global_config, plan = _load_project_plan(args.workflow, args.config)
    if args.project_command == "validate":
        print(
            f"Valid project workflow: {workflow.source_path} "
            f"({len(plan.stages)} stage(s))"
        )
        return 0
    if args.project_command == "plan":
        print(render_plan(plan, format="json" if args.json else "human"))
        return 0
    if args.project_command != "run":
        raise TxSuiteError(f"Unknown project command: {args.project_command}")

    selected = _selected_project_stages(args.stages)
    # Expand dependencies before a dry run returns and before an actual run writes.
    effective_stage_ids, auto_added = _effective_project_stages(
        plan,
        selected_stage_ids=selected,
        from_stage=args.from_stage,
        to_stage=args.to_stage,
    )
    if args.dry_run:
        print(_render_project_run_plan(plan, effective_stage_ids, auto_added))
        return 0

    runs_dir = _project_runs_dir(workflow)
    if args.resume:
        if args.run_id is None:
            raise TxSuiteError("--resume requires --run-id to identify an existing run")
        run_dir = runs_dir / args.run_id
        if (
            run_dir.parent != runs_dir
            or run_dir.resolve(strict=False).parent
            != runs_dir.resolve(strict=False)
        ):
            raise TxSuiteError(f"Unsafe run ID: {args.run_id!r}")
        if not run_dir.is_dir():
            raise TxSuiteError(f"Cannot resume missing run bundle: {run_dir}")
        bundle, manifest = _load_project_bundle(run_dir)
        if bundle.run_id != args.run_id:
            raise TxSuiteError(
                f"Cannot resume run {args.run_id!r}: manifest identifies "
                f"run {bundle.run_id!r}"
            )
        _require_compatible_resume(
            bundle, manifest, workflow, global_config, plan
        )
    else:
        bundle = RunBundle.create(
            runs_dir,
            run_id=args.run_id,
            workflow=workflow.to_dict(),
            resolved_config=global_config,
            command_plan=plan,
        )

    try:
        ProjectExecutor(bundle).execute(
            plan.to_dict(),
            resume=args.resume,
            config_hash=hash_payload(global_config),
            selected_stage_ids=effective_stage_ids,
            raise_on_failure=True,
        )
    except StageExecutionError as exc:
        raise TxSuiteError(f"{exc}; run directory: {exc.run_dir}") from exc
    print(bundle.run_dir)
    print(f"Selected stages: {', '.join(effective_stage_ids)}")
    if auto_added:
        print(f"Auto-added prerequisites: {', '.join(auto_added)}")
    return 0


def run(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "tools":
            for tool in select_tools(args.modality, args.stage):
                print("\t".join(tool))
            return 0
        if args.command == "config" and args.config_command == "show":
            print(json.dumps(load_config(args.file), indent=2, sort_keys=True))
            return 0
        if args.command == "config" and args.config_command == "init":
            if args.file.exists() and not args.force:
                raise ConfigError(f"{args.file} exists; use --force to replace it")
            args.file.parent.mkdir(parents=True, exist_ok=True)
            args.file.write_text(DEFAULT_TOML, encoding="utf-8")
            print(args.file)
            return 0
        if args.command == "project":
            return _run_project_command(args)
        if args.command == "workflow" and args.workflow_command == "bulk":
            config = load_config(args.config)
            validate_bulk_samplesheet(args.input)
            command = bulk_workflow_command(
                config,
                samplesheet=args.input,
                outdir=args.outdir,
                params_file=args.params_file,
                nextflow_config=args.nextflow_config,
                resume=args.resume,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            run_dir = args.run_dir or args.outdir / ".txsuite"
            run_command(
                command,
                run_dir=run_dir,
                task="workflow.bulk",
                backend=config["pipelines"]["bulk"]["name"],
                inputs={"samplesheet": str(args.input.resolve())},
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "directory",
                        "label": "nf-core/rnaseq results",
                        "path": str(args.outdir.resolve()),
                    }
                ],
            )
            return 0
        if args.command == "workflow" and args.workflow_command == "single-cell":
            config = load_config(args.config)
            validate_single_cell_samplesheet(args.input)
            command = single_cell_workflow_command(
                config,
                samplesheet=args.input,
                outdir=args.outdir,
                aligner=args.aligner,
                protocol=args.protocol,
                params_file=args.params_file,
                nextflow_config=args.nextflow_config,
                resume=args.resume,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="workflow.single-cell",
                backend=config["pipelines"]["single_cell"]["name"],
                inputs={
                    "samplesheet": str(args.input.resolve()),
                    "aligner": args.aligner,
                    "protocol": args.protocol,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "directory",
                        "label": "nf-core/scrnaseq results",
                        "path": str(args.outdir.resolve()),
                    }
                ],
            )
            return 0
        if args.command == "workflow" and args.workflow_command == "pseudobulk-de":
            config = load_config(args.config)
            command = pseudobulk_workflow_command(
                config,
                counts_layer=args.counts_layer,
                h5ad=args.input,
                outdir=args.outdir,
                sample_column=args.sample_column,
                design=args.design,
                reference=args.reference,
                test=args.test,
                group_column=args.group_column,
                group_value=args.group_value,
                covariates=tuple(args.covariate),
                single_cell_image=args.image,
                bulk_image=args.bulk_image,
                padj=args.padj,
                lfc=args.lfc,
                top_genes=args.top_genes,
                nextflow_config=args.nextflow_config,
                resume=args.resume,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="workflow.pseudobulk-de",
                backend="Nextflow DSL2",
                inputs={
                    "h5ad": str(args.input.resolve()),
                    "sample_column": args.sample_column,
                    "design": args.design,
                    "group": [args.group_column, args.group_value],
                    "covariates": args.covariate,
                    "contrast": [args.test, args.reference],
                    "padj": args.padj,
                    "abs_log2fc": args.lfc,
                    "top_genes": args.top_genes,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "directory",
                        "label": "pseudobulk tables",
                        "path": str((args.outdir / "pseudobulk").resolve()),
                    },
                    {
                        "kind": "directory",
                        "label": "DESeq2 results",
                        "path": str((args.outdir / "de" / f"{args.test}_vs_{args.reference}").resolve()),
                    },
                ],
            )
            return 0
        if args.command == "workflow" and args.workflow_command == "single-cell-cellranger":
            config = load_config(args.config)
            command = cellranger_workflow_command(
                config,
                genome_name=args.genome_name,
                fastqs=args.fastqs,
                sample=args.sample,
                outdir=args.outdir,
                fasta=args.fasta,
                gtf=args.gtf,
                reference=args.reference,
                cellranger_image=args.cellranger_image,
                threads=args.threads,
                memory_gb=args.memory_gb,
                create_bam=not args.no_bam,
                nextflow_config=args.nextflow_config,
                resume=args.resume,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="workflow.single-cell-cellranger",
                backend="Nextflow DSL2 + user-installed Cell Ranger",
                inputs={
                    "genome_name": args.genome_name,
                    "fastqs": str(args.fastqs.resolve()),
                    "sample": args.sample,
                    "fasta": str(args.fasta.resolve()) if args.fasta else None,
                    "gtf": str(args.gtf.resolve()) if args.gtf else None,
                    "reference": str(args.reference.resolve()) if args.reference else None,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "directory",
                        "label": "Cell Ranger counts",
                        "path": str((args.outdir / "counts" / args.sample / "outs").resolve()),
                    },
                    {
                        "kind": "directory",
                        "label": "filtered feature-barcode matrix",
                        "path": str(
                            (
                                args.outdir
                                / "counts"
                                / args.sample
                                / "outs"
                                / "filtered_feature_bc_matrix"
                            ).resolve()
                        ),
                    },
                ],
            )
            return 0
        if args.command == "workflow" and args.workflow_command == "spatial":
            config = load_config(args.config)
            command = spaceranger_command(
                run_id=args.id,
                transcriptome=args.transcriptome,
                fastqs=args.fastqs,
                image=args.image,
                sample=args.sample,
                slide=args.slide,
                area=args.area,
                unknown_slide=args.unknown_slide,
                create_bam=not args.no_bam,
                cores=args.cores,
                memory=args.memory,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            result_dir = args.outdir / args.id / "outs"
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / args.id / ".txsuite",
                task="workflow.spatial",
                backend=f"{config['pipelines']['spatial']['name']} {config['pipelines']['spatial']['release']}",
                inputs={
                    "transcriptome": str(args.transcriptome.resolve()),
                    "fastqs": str(args.fastqs.resolve()),
                    "image": str(args.image.resolve()),
                    "sample": args.sample,
                    "slide": args.slide,
                    "area": args.area,
                },
                outputs={"outdir": str(result_dir.resolve())},
                artifacts=[
                    {
                        "kind": "directory",
                        "label": "Space Ranger outputs",
                        "path": str(result_dir.resolve()),
                    },
                    {
                        "kind": "html",
                        "label": "Space Ranger web summary",
                        "path": str((result_dir / "web_summary.html").resolve()),
                    },
                ],
                cwd=args.outdir.resolve(),
            )
            return 0
        if args.command == "workflow" and args.workflow_command == "spatial-open":
            command = spacemake_command(args.project_root, args.cores)
            if args.dry_run:
                print(format_command(command))
                return 0
            run_command(
                command,
                run_dir=args.run_dir or args.project_root / ".txsuite",
                task="workflow.spatial-open",
                backend="Spacemake (experimental)",
                inputs={"project_root": str(args.project_root.resolve())},
                outputs={"project_root": str(args.project_root.resolve())},
                artifacts=[
                    {
                        "kind": "directory",
                        "label": "Spacemake project",
                        "path": str(args.project_root.resolve()),
                    }
                ],
                cwd=args.project_root.resolve(),
            )
            return 0
        if args.command == "bulk" and args.bulk_command == "de":
            config = load_config(args.config)
            image = args.image or config["images"]["bulk_r"]
            command = differential_expression_command(
                method=args.method,
                image=image,
                counts=args.counts,
                metadata=args.metadata,
                design=args.design,
                reference=args.reference,
                test=args.test,
                formula=args.formula,
                coefficient=args.coefficient,
                outdir=args.outdir,
                covariates=tuple(args.covariate),
                padj=args.padj,
                lfc=args.lfc,
                top_genes=args.top_genes,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            backend = {
                "deseq2": "DESeq2",
                "edger": "edgeR",
                "limma": "limma-voom",
            }[args.method]
            result_name = f"{args.method}-results.tsv"
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task=f"bulk.de.{args.method}",
                backend=backend,
                inputs={
                    "counts": str(args.counts.resolve()),
                    "metadata": str(args.metadata.resolve()),
                    "method": args.method,
                    "design": args.design or args.formula,
                    "contrast": (
                        [args.test, args.reference] if args.design else args.coefficient
                    ),
                    "covariates": args.covariate,
                    "padj": args.padj,
                    "abs_log2fc": args.lfc,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "table",
                        "label": f"{backend} results",
                        "path": str((args.outdir / result_name).resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "normalized counts",
                        "path": str((args.outdir / "normalized-counts.tsv").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "significant genes",
                        "path": str((args.outdir / "significant-genes.tsv").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "sample QC",
                        "path": str((args.outdir / "sample-qc.tsv").resolve()),
                    },
                    {
                        "kind": "figure",
                        "label": "sample ordination",
                        "path": str(
                            (
                                args.outdir
                                / ("pca.pdf" if args.method == "deseq2" else "mds.pdf")
                            ).resolve()
                        ),
                    },
                    {
                        "kind": "figure",
                        "label": "volcano plot",
                        "path": str((args.outdir / "volcano.pdf").resolve()),
                    },
                    {
                        "kind": "text",
                        "label": "R session info",
                        "path": str((args.outdir / "session-info.txt").resolve()),
                    },
                ],
            )
            return 0
        if args.command == "bulk" and args.bulk_command == "enrich":
            config = load_config(args.config)
            image = args.image or config["images"]["bulk_r"]
            command = enrichment_command(
                image=image,
                de_results=args.de,
                genesets=args.genesets,
                mode=args.mode,
                outdir=args.outdir,
                padj=args.padj,
                lfc=args.lfc,
                min_size=args.min_size,
                max_size=args.max_size,
                adjust=args.adjust,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            result_name = (
                "ora-results.tsv" if args.mode == "ora" else "gsea-results.tsv"
            )
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task=f"bulk.enrich.{args.mode}",
                backend="clusterProfiler",
                inputs={
                    "de_results": str(args.de.resolve()),
                    "genesets": str(args.genesets.resolve()),
                    "mode": args.mode,
                    "padj": args.padj,
                    "abs_log2fc": args.lfc,
                    "gene_set_size": [args.min_size, args.max_size],
                    "adjust": args.adjust,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "table",
                        "label": f"{args.mode.upper()} results",
                        "path": str((args.outdir / result_name).resolve()),
                    },
                    {
                        "kind": "figure",
                        "label": "top enrichment results",
                        "path": str((args.outdir / "enrichment-top.pdf").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "enrichment summary",
                        "path": str((args.outdir / "enrichment-summary.tsv").resolve()),
                    },
                ],
            )
            return 0
        if args.command == "single-cell" and args.single_cell_command == "analyze":
            config = load_config(args.config)
            image = args.image or config["images"]["single_cell_python"]
            command = analysis_command(
                image=image,
                input_path=args.input,
                outdir=args.outdir,
                min_genes=args.min_genes,
                min_cells=args.min_cells,
                max_mito_pct=args.max_mito_pct,
                resolution=args.resolution,
                metadata=args.metadata,
                barcode_column=args.barcode_column,
                batch_column=args.batch_column,
                integration=args.integration,
                counts_layer=args.counts_layer,
                target_sum=args.target_sum,
                n_hvg=args.n_hvg,
                hvg_flavor=args.hvg_flavor,
                n_pcs=args.n_pcs,
                n_neighbors=args.n_neighbors,
                umap_min_dist=args.umap_min_dist,
                marker_method=args.marker_method,
                stop_after=args.stop_after,
                skip_umap=args.skip_umap,
                skip_markers=args.skip_markers,
                doublets=args.doublets,
                doublet_batch_column=args.doublet_batch_column,
                expected_doublet_rate=args.expected_doublet_rate,
                doublet_threshold=args.doublet_threshold,
                top_markers=args.top_markers,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="single-cell.analyze",
                backend="Scanpy",
                inputs={
                    "data": str(args.input.resolve()),
                    "metadata": str(args.metadata.resolve()) if args.metadata else None,
                    "barcode_column": args.barcode_column,
                    "batch_column": args.batch_column,
                    "integration": args.integration,
                    "counts_layer": args.counts_layer,
                    "target_sum": args.target_sum,
                    "n_hvg": args.n_hvg,
                    "hvg_flavor": args.hvg_flavor,
                    "n_pcs": args.n_pcs,
                    "n_neighbors": args.n_neighbors,
                    "umap_min_dist": args.umap_min_dist,
                    "marker_method": args.marker_method,
                    "stop_after": args.stop_after,
                    "skip_umap": args.skip_umap,
                    "skip_markers": args.skip_markers,
                    "min_genes": args.min_genes,
                    "min_cells": args.min_cells,
                    "max_mito_pct": args.max_mito_pct,
                    "resolution": args.resolution,
                    "doublets": args.doublets,
                    "doublet_batch_column": args.doublet_batch_column,
                    "expected_doublet_rate": args.expected_doublet_rate,
                    "doublet_threshold": args.doublet_threshold,
                    "top_markers": args.top_markers,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "anndata",
                        "label": "analyzed AnnData",
                        "path": str((args.outdir / "analysis.h5ad").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "cell QC metrics",
                        "path": str((args.outdir / "cell-qc.tsv").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "Leiden clusters and UMAP",
                        "path": str((args.outdir / "clusters.tsv").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "Leiden marker genes",
                        "path": str((args.outdir / "marker-genes.tsv").resolve()),
                    },
                ],
            )
            return 0
        if (
            args.command == "single-cell"
            and args.single_cell_command == "pseudobulk-batch"
        ):
            config = load_config(args.config)
            if args.direct:
                commands = run_pseudobulk_manifest(
                    counts_layer=args.counts_layer,
                    manifest=args.manifest,
                    h5ad=args.input,
                    outdir=args.outdir,
                    sample_column=args.sample_column,
                    single_cell_image=args.image or config["images"]["single_cell_python"],
                    bulk_image=args.bulk_image or config["images"]["bulk_r"],
                    resume=args.resume,
                    dry_run=args.dry_run,
                )
                for command in commands:
                    print(format_command(command))
                return 0
            command = pseudobulk_manifest_workflow_command(
                config,
                counts_layer=args.counts_layer,
                manifest=args.manifest,
                h5ad=args.input,
                outdir=args.outdir,
                sample_column=args.sample_column,
                single_cell_image=args.image,
                bulk_image=args.bulk_image,
                nextflow_config=args.nextflow_config,
                resume=args.resume,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="single-cell.pseudobulk-batch",
                backend="Nextflow DSL2",
                inputs={
                    "h5ad": str(args.input.resolve()),
                    "manifest": str(args.manifest.resolve()),
                    "sample_column": args.sample_column,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {"kind": "table", "label": "comparison index", "path": str((args.outdir / "comparison-index.tsv").resolve())},
                    {"kind": "table", "label": "combined DE results", "path": str((args.outdir / "combined-results.tsv").resolve())},
                ],
            )
            return 0
        if (
            args.command == "single-cell"
            and args.single_cell_command == "pseudobulk-de"
        ):
            config = load_config(args.config)
            image = args.image or config["images"]["single_cell_python"]
            bulk_image = args.bulk_image or config["images"]["bulk_r"]
            counts = args.outdir / "pseudobulk-counts.tsv"
            metadata = args.outdir / "pseudobulk-metadata.tsv"
            aggregate = pseudobulk_command(
                counts_layer=args.counts_layer,
                image=image,
                h5ad=args.input,
                outdir=args.outdir,
                sample_column=args.sample_column,
                design=args.design,
                group_column=args.group_column,
                group_value=args.group_value,
                covariates=tuple(args.covariate),
                reference=args.reference,
                test=args.test,
            )
            if args.dry_run:
                differential = deseq2_command(
                    image=bulk_image,
                    counts=counts,
                    metadata=metadata,
                    design=args.design,
                    reference=args.reference,
                    test=args.test,
                    outdir=args.outdir,
                    covariates=tuple(args.covariate),
                    check_inputs=False,
                )
                print(format_command(aggregate))
                print(format_command(differential))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            base_run_dir = args.run_dir or args.outdir / ".txsuite"
            run_command(
                aggregate,
                run_dir=base_run_dir / "pseudobulk",
                task="single-cell.pseudobulk",
                backend="Scanpy",
                inputs={
                    "h5ad": str(args.input.resolve()),
                    "sample_column": args.sample_column,
                    "design": args.design,
                    "group": [args.group_column, args.group_value],
                    "covariates": args.covariate,
                },
                outputs={
                    "counts": str(counts.resolve()),
                    "metadata": str(metadata.resolve()),
                },
                artifacts=[
                    {
                        "kind": "table",
                        "label": "pseudobulk counts",
                        "path": str(counts.resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "pseudobulk metadata",
                        "path": str(metadata.resolve()),
                    },
                ],
            )
            differential = deseq2_command(
                image=bulk_image,
                counts=counts,
                metadata=metadata,
                design=args.design,
                reference=args.reference,
                test=args.test,
                outdir=args.outdir,
                covariates=tuple(args.covariate),
            )
            run_command(
                differential,
                run_dir=base_run_dir / "deseq2",
                task="single-cell.pseudobulk-de",
                backend="DESeq2",
                inputs={
                    "counts": str(counts.resolve()),
                    "metadata": str(metadata.resolve()),
                    "design": args.design,
                    "group": [args.group_column, args.group_value],
                    "covariates": args.covariate,
                    "contrast": [args.test, args.reference],
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "table",
                        "label": "DESeq2 results",
                        "path": str((args.outdir / "deseq2-results.tsv").resolve()),
                    }
                ],
            )
            return 0
        if args.command == "spatial" and args.spatial_command == "analyze":
            config = load_config(args.config)
            image = args.image or config["images"]["spatial_python"]
            command = spatial_analysis_command(
                image=image,
                input_path=args.input,
                outdir=args.outdir,
                dataset_id=args.dataset_id,
                min_counts=args.min_counts,
                min_spots=args.min_spots,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="spatial.analyze",
                backend="SpatialData/Squidpy",
                inputs={
                    "data": str(args.input.resolve()),
                    "dataset_id": args.dataset_id,
                    "min_counts": args.min_counts,
                    "min_spots": args.min_spots,
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "spatialdata",
                        "label": "SpatialData Zarr",
                        "path": str((args.outdir / "spatialdata.zarr").resolve()),
                    },
                    {
                        "kind": "anndata",
                        "label": "spatial analysis",
                        "path": str((args.outdir / "analysis.h5ad").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "spot QC metrics",
                        "path": str((args.outdir / "spot-qc.tsv").resolve()),
                    },
                ],
            )
            return 0
        if args.command == "reference" and args.reference_command == "cache":
            print(cache_reference(args.source, args.sha256, args.name, args.root))
            return 0
        if args.command == "env" and args.env_command == "list":
            for name, image in load_config(args.config)["images"].items():
                print(f"{name.replace('_', '-')}\t{image}")
            return 0
        if args.command == "env" and args.env_command == "doctor":
            config = load_config(args.config)
            profile = config["execution"]["profile"]
            required = ("nextflow", profile)
            optional = tuple(
                executable
                for executable in (
                    "docker",
                    "cellranger",
                    "spaceranger",
                    "spacemake",
                )
                if executable not in required
            )
            missing = False
            for executable in (*required, *optional):
                location = shutil.which(executable)
                label = "FOUND" if location else "MISSING"
                requirement = "required" if executable in required else "optional"
                print(f"{label}\t{requirement}\t{executable}\t{location or '-'}")
                missing |= executable in required and location is None
            return int(missing)
        if args.command == "env" and args.env_command == "verify-images":
            config = load_config(args.config)
            mutable = False
            for name, image in config["images"].items():
                locked = image_is_locked(image)
                print(f"{'LOCKED' if locked else 'MUTABLE'}\t{name}\t{image}")
                mutable |= not locked
            return int(mutable)
        if args.command == "env" and args.env_command == "build":
            if args.environment == "cellranger":
                if not args.tag:
                    raise TxSuiteError(
                        "--tag is required for 'cellranger'; TxSuite has no default "
                        "image because it never publishes this licensed software"
                    )
                if not args.source_tarball:
                    raise TxSuiteError(
                        "--source-tarball is required for 'cellranger': point it at "
                        "your own EULA-accepted cellranger-*.tar.gz download"
                    )
                tag = args.tag
                if args.dry_run:
                    print(
                        f"docker build --tag {shlex.quote(tag)} "
                        f"--build-arg CELLRANGER_TARBALL={shlex.quote(args.source_tarball.name)} "
                        "<bundled-cellranger-context>"
                    )
                    return 0
                run_dir = args.run_dir or Path(".txsuite/build-cellranger")
                build_cellranger_image(
                    tag, source_tarball=args.source_tarball, run_dir=run_dir
                )
                return 0
            config = load_config(args.config)
            key = args.environment.replace("-", "_")
            tag = args.tag or config["images"][key]
            if args.dry_run:
                print(
                    f"docker build --tag {shlex.quote(tag)} "
                    f"<bundled-{args.environment}-context>"
                )
                return 0
            run_dir = args.run_dir or Path(f".txsuite/build-{args.environment}")
            if args.environment == "bulk-r":
                build_bulk_r_image(tag, run_dir=run_dir)
            elif args.environment == "single-cell-python":
                build_single_cell_image(tag, run_dir=run_dir)
            else:
                build_spatial_image(tag, run_dir=run_dir)
            return 0
    except (
        ConfigError,
        WorkflowConfigError,
        ProvenanceError,
        TxSuiteError,
        OSError,
        ValueError,
        KeyError,
    ) as exc:
        print(f"error: {exc}")
        return 2
    return 0


def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()
