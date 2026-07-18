from __future__ import annotations

import argparse
import json
import shlex
import shutil
from pathlib import Path

from txsuite.bulk import (
    build_bulk_r_image,
    deseq2_command,
    validate_samplesheet as validate_bulk_samplesheet,
    workflow_command as bulk_workflow_command,
)
from txsuite.catalog import select_tools
from txsuite.config import ConfigError, DEFAULT_TOML, load_config
from txsuite.hardening import cache_reference, image_is_locked
from txsuite.runtime import TxSuiteError, format_command, run_command
from txsuite.single_cell import (
    analysis_command,
    build_single_cell_image,
    pseudobulk_command,
    validate_samplesheet as validate_single_cell_samplesheet,
    workflow_command as single_cell_workflow_command,
)
from txsuite.spatial import (
    analysis_command as spatial_analysis_command,
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
    de = bulk_commands.add_parser("de", help="run DESeq2 in the bulk-r image")
    de.add_argument("--counts", type=Path, required=True)
    de.add_argument("--metadata", type=Path, required=True)
    de.add_argument("--design", required=True)
    de.add_argument("--reference", required=True)
    de.add_argument("--test", required=True)
    de.add_argument("--outdir", type=Path, required=True)
    de.add_argument("--image")
    de.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    de.add_argument("--run-dir", type=Path)
    de.add_argument("--dry-run", action="store_true")

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
    analyze.add_argument("--image")
    analyze.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    analyze.add_argument("--run-dir", type=Path)
    analyze.add_argument("--dry-run", action="store_true")

    pseudobulk = single_cell_commands.add_parser(
        "pseudobulk-de", help="aggregate cells and run the bulk DESeq2 backend"
    )
    pseudobulk.add_argument("--input", type=Path, required=True)
    pseudobulk.add_argument("--sample-column", required=True)
    pseudobulk.add_argument("--design", required=True)
    pseudobulk.add_argument("--reference", required=True)
    pseudobulk.add_argument("--test", required=True)
    pseudobulk.add_argument("--outdir", type=Path, required=True)
    pseudobulk.add_argument("--image")
    pseudobulk.add_argument("--bulk-image")
    pseudobulk.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    pseudobulk.add_argument("--run-dir", type=Path)
    pseudobulk.add_argument("--dry-run", action="store_true")

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
    doctor = env_commands.add_parser("doctor", help="check workflow prerequisites")
    doctor.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    verify_images = env_commands.add_parser(
        "verify-images", help="fail when configured images use mutable tags"
    )
    verify_images.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    build = env_commands.add_parser("build", help="build a TxSuite-owned image")
    build.add_argument(
        "environment", choices=("bulk-r", "single-cell-python", "spatial-python")
    )
    build.add_argument("--tag")
    build.add_argument("--config", type=Path, default=Path("txsuite.toml"))
    build.add_argument("--run-dir", type=Path)
    build.add_argument("--dry-run", action="store_true")
    return parser


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
            command = deseq2_command(
                image=image,
                counts=args.counts,
                metadata=args.metadata,
                design=args.design,
                reference=args.reference,
                test=args.test,
                outdir=args.outdir,
            )
            if args.dry_run:
                print(format_command(command))
                return 0
            args.outdir.mkdir(parents=True, exist_ok=True)
            run_command(
                command,
                run_dir=args.run_dir or args.outdir / ".txsuite",
                task="bulk.de",
                backend="DESeq2",
                inputs={
                    "counts": str(args.counts.resolve()),
                    "metadata": str(args.metadata.resolve()),
                    "design": args.design,
                    "contrast": [args.test, args.reference],
                },
                outputs={"outdir": str(args.outdir.resolve())},
                artifacts=[
                    {
                        "kind": "table",
                        "label": "DESeq2 results",
                        "path": str((args.outdir / "deseq2-results.tsv").resolve()),
                    },
                    {
                        "kind": "table",
                        "label": "normalized counts",
                        "path": str((args.outdir / "normalized-counts.tsv").resolve()),
                    },
                    {
                        "kind": "text",
                        "label": "R session info",
                        "path": str((args.outdir / "session-info.txt").resolve()),
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
                    "min_genes": args.min_genes,
                    "min_cells": args.min_cells,
                    "max_mito_pct": args.max_mito_pct,
                    "resolution": args.resolution,
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
                image=image,
                h5ad=args.input,
                outdir=args.outdir,
                sample_column=args.sample_column,
                design=args.design,
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
    except (ConfigError, TxSuiteError) as exc:
        print(f"error: {exc}")
        return 2
    return 0


def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()
