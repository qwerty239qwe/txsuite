"""Stage adapters backed by :mod:`txsuite.single_cell` command builders."""

from __future__ import annotations

import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from txsuite.single_cell import (
    alevin_workflow_command,
    analysis_command,
    pseudobulk_command as _pseudobulk_command,
    pseudobulk_workflow_command,
    workflow_command,
)

from . import (
    command_context,
    configured_image,
    input_path,
    optional_input_path,
)


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
        simpleaf_index=_optional_path(params.get("simpleaf_index")),
        txp2gene=_optional_path(params.get("txp2gene")),
        resume=resume,
    )


def alevin_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the native simpleaf/alevin-fry DAG for a registry context."""

    config, inputs, params, outdir, resume, check_inputs = command_context(context)
    return alevin_workflow_command(
        dict(config),
        samplesheet=input_path(inputs, "samplesheet"),
        outdir=outdir,
        fasta=_optional_path(params.get("fasta")),
        gtf=_optional_path(params.get("gtf")),
        simpleaf_index=_optional_path(params.get("simpleaf_index")),
        chemistry=params.get("chemistry", "10xv3"),
        resolution=params.get("resolution", "cr-like"),
        whitelist=_optional_path(params.get("whitelist")),
        rlen=params.get("rlen", 91),
        salmon_image=params.get("image"),
        single_cell_image=params.get("single_cell_image"),
        nextflow_config=_optional_path(params.get("nextflow_config")),
        resume=resume,
        check_inputs=check_inputs,
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
        metadata=optional_input_path(inputs, "metadata"),
        barcode_column=params.get("barcode_column", "barcode"),
        batch_column=params.get("batch_column"),
        integration=params.get("integration", "none"),
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


build_alevin_command = alevin_command
build_scrnaseq_command = scrnaseq_command
build_scanpy_command = scanpy_command
build_pseudobulk_command = pseudobulk_command
build_pseudobulk_de_command = pseudobulk_de_command


__all__ = [
    "alevin_command",
    "build_alevin_command",
    "build_pseudobulk_command",
    "build_pseudobulk_de_command",
    "build_scanpy_command",
    "build_scrnaseq_command",
    "pseudobulk_command",
    "pseudobulk_de_command",
    "scanpy_command",
    "scrnaseq_command",
]
