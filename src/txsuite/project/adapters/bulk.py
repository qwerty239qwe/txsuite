"""Stage adapters backed by :mod:`txsuite.bulk` command builders."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from txsuite.bulk import (
    differential_expression_command,
    enrichment_command,
    salmon_workflow_command,
    workflow_command,
)

from . import command_context, configured_image, input_path


def _optional_path(value: Any) -> Path | None:
    return None if value is None else Path(value)


def bulk_rnaseq_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the pinned nf-core/rnaseq command for a registry context."""

    config, inputs, params, outdir, resume, _ = command_context(context)
    return workflow_command(
        dict(config),
        samplesheet=input_path(inputs, "samplesheet"),
        outdir=outdir,
        params_file=_optional_path(params.get("params_file")),
        nextflow_config=_optional_path(params.get("nextflow_config")),
        pseudo_aligner=params.get("pseudo_aligner"),
        skip_alignment=params.get("skip_alignment", False),
        salmon_index=_optional_path(params.get("salmon_index")),
        resume=resume,
    )


def bulk_salmon_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build the native salmon quantification DAG for a registry context."""

    config, inputs, params, outdir, resume, check_inputs = command_context(context)
    return salmon_workflow_command(
        dict(config),
        samplesheet=input_path(inputs, "samplesheet"),
        outdir=outdir,
        fasta=_optional_path(params.get("fasta")),
        gtf=_optional_path(params.get("gtf")),
        salmon_index=_optional_path(params.get("salmon_index")),
        tx2gene=_optional_path(params.get("tx2gene")),
        salmon_image=params.get("image"),
        libtype=params.get("libtype", "A"),
        kmer_len=params.get("kmer_len", 31),
        gencode=params.get("gencode", False),
        nextflow_config=_optional_path(params.get("nextflow_config")),
        resume=resume,
        check_inputs=check_inputs,
    )


def bulk_de_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build a bulk differential-expression container command."""

    config, inputs, params, outdir, _, check_inputs = command_context(context)
    return differential_expression_command(
        method=params.get("method", "deseq2"),
        image=configured_image(config, params, "bulk_r"),
        counts=input_path(inputs, "counts"),
        metadata=input_path(inputs, "metadata"),
        design=params["design"],
        reference=params["reference"],
        test=params["test"],
        outdir=outdir,
        covariates=tuple(params.get("covariates", ())),
        padj=params.get("padj", 0.05),
        lfc=params.get("lfc", 1.0),
        top_genes=params.get("top_genes", 50),
        contrasts=params.get("contrasts", "single"),
        check_inputs=check_inputs,
    )


def bulk_enrichment_command(context: Mapping[str, Any] | object) -> list[str]:
    """Build a bulk ORA/GSEA container command."""

    config, inputs, params, outdir, _, check_inputs = command_context(context)
    return enrichment_command(
        image=configured_image(config, params, "bulk_r"),
        de_results=input_path(inputs, "de_results"),
        genesets=input_path(inputs, "genesets"),
        mode=params.get("mode", "ora"),
        outdir=outdir,
        padj=params.get("padj", 0.05),
        lfc=params.get("lfc", 1.0),
        min_size=params.get("min_size", 10),
        max_size=params.get("max_size", 500),
        adjust=params.get("adjust", "BH"),
        check_inputs=check_inputs,
    )


# Verbose aliases are useful to callers that distinguish factories from commands.
build_bulk_rnaseq_command = bulk_rnaseq_command
build_bulk_salmon_command = bulk_salmon_command
build_bulk_de_command = bulk_de_command
build_bulk_enrichment_command = bulk_enrichment_command
rnaseq_command = bulk_rnaseq_command
de_command = bulk_de_command
enrichment_stage_command = bulk_enrichment_command


__all__ = [
    "build_bulk_de_command",
    "build_bulk_enrichment_command",
    "build_bulk_rnaseq_command",
    "build_bulk_salmon_command",
    "bulk_de_command",
    "bulk_enrichment_command",
    "bulk_rnaseq_command",
    "bulk_salmon_command",
    "de_command",
    "enrichment_stage_command",
    "rnaseq_command",
]
