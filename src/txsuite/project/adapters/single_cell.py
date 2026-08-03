"""Stage adapters backed by :mod:`txsuite.single_cell` command builders."""

from __future__ import annotations

import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from txsuite.single_cell import (
    analysis_command,
    pseudobulk_command as _pseudobulk_command,
    pseudobulk_workflow_command,
    workflow_command,
)

from . import command_context, configured_image, input_path


def _optional_path(value: Any) -> Path | None:
    return None if value is None else Path(value)


def scrnaseq_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the pinned nf-core/scrnaseq command for a registry context."""

    config, inputs, params, outdir, resume, _ = command_context(context)
    return workflow_command(
        dict(config),
        samplesheet=input_path(inputs, "samplesheet"),
        outdir=outdir,
        aligner=params.get("aligner", "simpleaf"),
        protocol=params.get("protocol"),
        params_file=_optional_path(params.get("params_file")),
        nextflow_config=_optional_path(params.get("nextflow_config")),
        resume=resume,
    )


def scanpy_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the TxSuite Scanpy analysis command."""

    config, inputs, params, outdir, _, check_inputs = command_context(context)
    planned_input = input_path(inputs, "input")
    symbolic = bool(re.fullmatch(r"\$\{[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\}", str(planned_input)))
    return analysis_command(
        image=configured_image(config, params, "single_cell_python"),
        input_path=planned_input,
        outdir=outdir,
        min_genes=params.get("min_genes", 200),
        min_cells=params.get("min_cells", 3),
        max_mito_pct=params.get("max_mito_pct", 20.0),
        resolution=params.get("resolution", 1.0),
        check_inputs=check_inputs,
        container_target="/input/data.h5ad" if symbolic else None,
    )


def pseudobulk_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the single-cell-to-pseudobulk aggregation command."""

    config, inputs, params, outdir, _, check_inputs = command_context(context)
    return _pseudobulk_command(
        image=configured_image(config, params, "single_cell_python"),
        h5ad=input_path(inputs, "h5ad"),
        outdir=outdir,
        sample_column=params["sample_column"],
        design=params["design"],
        check_inputs=check_inputs,
    )


def pseudobulk_de_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the native resumable pseudobulk differential-expression DAG."""

    config, inputs, params, outdir, resume, check_inputs = command_context(context)
    return pseudobulk_workflow_command(
        dict(config),
        h5ad=input_path(inputs, "h5ad"),
        outdir=outdir,
        sample_column=params["sample_column"],
        design=params["design"],
        reference=params["reference"],
        test=params["test"],
        single_cell_image=params.get("single_cell_image"),
        bulk_image=params.get("bulk_image"),
        padj=params.get("padj", 0.05),
        lfc=params.get("lfc", 1.0),
        top_genes=params.get("top_genes", 50),
        nextflow_config=_optional_path(params.get("nextflow_config")),
        resume=resume,
        check_inputs=check_inputs,
    )


build_scrnaseq_command = scrnaseq_command
build_scanpy_command = scanpy_command
build_pseudobulk_command = pseudobulk_command
build_pseudobulk_de_command = pseudobulk_de_command


__all__ = [
    "build_pseudobulk_command",
    "build_pseudobulk_de_command",
    "build_scanpy_command",
    "build_scrnaseq_command",
    "pseudobulk_command",
    "pseudobulk_de_command",
    "scanpy_command",
    "scrnaseq_command",
]
