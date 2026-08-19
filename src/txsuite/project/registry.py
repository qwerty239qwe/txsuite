"""Declarative registry of project-stage backends."""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

from txsuite.alignment import STRANDEDNESS
from txsuite.bulk import CONTRAST_MODES, PSEUDO_ALIGNERS, SALMON_LIBTYPES
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import (
    ALEVIN_CHEMISTRIES,
    ALEVIN_RESOLUTIONS,
    INTEGRATION_METHODS,
)
from txsuite.variants import ANNOTATION_TOOLS

from .adapters.bulk import (
    bulk_align_command,
    bulk_de_command,
    bulk_enrichment_command,
    bulk_rnaseq_command,
    bulk_rnavar_command,
    bulk_salmon_command,
    bulk_star_reference_command,
)
from .adapters.single_cell import (
    alevin_command,
    pseudobulk_command,
    pseudobulk_de_command,
    scanpy_command,
    scrnaseq_command,
)


StageCommandFactory = Callable[[Mapping[str, Any] | object], list[str]]
ParameterValidator = Callable[[Any], Any]


@dataclass(frozen=True, slots=True)
class OutputPolicy:
    """Filesystem verification policy for one declared stage output."""

    kind: str
    non_empty: bool = False

    def __post_init__(self) -> None:
        if self.kind not in {"file", "directory"}:
            raise ValueError("OutputPolicy kind must be 'file' or 'directory'")
        if not isinstance(self.non_empty, bool):
            raise TypeError("OutputPolicy non_empty must be a boolean")

    def to_dict(self) -> dict[str, Any]:
        return {"kind": self.kind, "non_empty": self.non_empty}


@dataclass(frozen=True)
class StageSpec:
    """Immutable declaration of one runnable TxSuite stage.

    ``id`` and ``uses`` intentionally name the same canonical key.  ``id`` is the
    concise registry spelling while ``uses`` mirrors the project-file field used by
    planners and runners.
    """

    id: str
    modality: str
    maturity: str
    inputs: Mapping[str, str]
    outputs: Mapping[str, str]
    defaults: Mapping[str, Any]
    validators: Mapping[str, ParameterValidator]
    command_factory: StageCommandFactory
    required_executables: tuple[str, ...] = ()
    required_images: tuple[str, ...] = ()
    supports_resume: bool = False
    uses: str | None = None
    output_policies: Mapping[str, OutputPolicy | Mapping[str, Any]] = field(
        default_factory=dict
    )
    optional_inputs: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        canonical = self.id if self.uses is None else self.uses
        if canonical != self.id:
            raise ValueError("StageSpec id and uses must contain the same canonical key")
        if not re.fullmatch(r"[a-z][a-z0-9-]*(?:\.[a-z][a-z0-9-]*)+", self.id):
            raise ValueError(f"Invalid stage id: {self.id!r}")
        if self.maturity not in {"ready", "selected", "experimental", "external"}:
            raise ValueError(f"Invalid stage maturity: {self.maturity!r}")
        if not callable(self.command_factory):
            raise TypeError("StageSpec command_factory must be callable")
        object.__setattr__(self, "uses", canonical)
        object.__setattr__(self, "inputs", _freeze_contract(self.inputs, "input"))
        optional_inputs = _freeze_contract(self.optional_inputs, "optional input")
        overlapping = sorted(set(optional_inputs) & set(self.inputs))
        if overlapping:
            raise ValueError(
                "Stage inputs cannot be both required and optional: "
                + ", ".join(overlapping)
            )
        object.__setattr__(self, "optional_inputs", optional_inputs)
        outputs = _freeze_contract(self.outputs, "output")
        object.__setattr__(self, "outputs", outputs)
        object.__setattr__(
            self,
            "output_policies",
            _freeze_output_policies(outputs, self.output_policies),
        )
        object.__setattr__(self, "defaults", _deep_freeze(self.defaults))
        object.__setattr__(self, "validators", MappingProxyType(dict(self.validators)))
        object.__setattr__(self, "required_executables", tuple(self.required_executables))
        object.__setattr__(self, "required_images", tuple(self.required_images))

    @property
    def accepted_inputs(self) -> Mapping[str, str]:
        """Every input this stage understands, required and optional alike."""

        return MappingProxyType({**self.inputs, **self.optional_inputs})

    def build_command(self, context: Mapping[str, Any] | object) -> list[str]:
        """Validate/default context parameters, then invoke the command factory."""

        if isinstance(context, Mapping):
            raw_params = context.get("params", {})
            normalized = validate_stage_parameters(self, raw_params)
            resolved = dict(context)
            resolved["params"] = normalized
            return self.command_factory(resolved)
        raw_params = getattr(context, "params", {})
        normalized = validate_stage_parameters(self, raw_params)
        values = {
            "global_config": getattr(
                context, "global_config", getattr(context, "config", None)
            ),
            "resolved_inputs": getattr(
                context, "resolved_inputs", getattr(context, "inputs", None)
            ),
            "params": normalized,
            "outdir": getattr(context, "outdir", None),
            "resume": getattr(context, "resume", False),
            "check_inputs": getattr(context, "check_inputs", True),
        }
        return self.command_factory(values)


def _freeze_contract(contract: Mapping[str, str], label: str) -> Mapping[str, str]:
    values = dict(contract)
    for name, artifact_type in values.items():
        if not isinstance(name, str) or not name:
            raise ValueError(f"Stage {label} names must be non-empty strings")
        if not isinstance(artifact_type, str) or not artifact_type:
            raise ValueError(f"Stage {label} artifact types must be non-empty strings")
        if any(character in artifact_type for character in "*?[]/\\"):
            raise ValueError(
                f"Stage {label} contract {name!r} must name an artifact type, not a glob"
            )
    return MappingProxyType(values)


def _freeze_output_policies(
    outputs: Mapping[str, str],
    policies: Mapping[str, OutputPolicy | Mapping[str, Any]],
) -> Mapping[str, OutputPolicy]:
    supplied = dict(policies)
    if not supplied:
        # Backward compatibility for third-party StageSpec construction. Built-in
        # stages declare every policy explicitly below.
        supplied = {
            name: OutputPolicy(
                kind="directory" if name == "results" else "file",
                non_empty=False,
            )
            for name in outputs
        }
    if set(supplied) != set(outputs):
        missing = sorted(set(outputs) - set(supplied))
        extra = sorted(set(supplied) - set(outputs))
        details: list[str] = []
        if missing:
            details.append(f"missing: {', '.join(missing)}")
        if extra:
            details.append(f"unknown: {', '.join(extra)}")
        raise ValueError("Stage output policies must exactly match outputs (" + "; ".join(details) + ")")
    normalized: dict[str, OutputPolicy] = {}
    for name, value in supplied.items():
        if isinstance(value, OutputPolicy):
            policy = value
        elif isinstance(value, Mapping):
            unknown = sorted(set(value) - {"kind", "non_empty"})
            if unknown:
                raise ValueError(
                    f"Output policy {name!r} has unknown key(s): {', '.join(unknown)}"
                )
            if "kind" not in value:
                raise ValueError(f"Output policy {name!r} requires kind")
            policy = OutputPolicy(
                kind=value["kind"], non_empty=value.get("non_empty", False)
            )
        else:
            raise TypeError(f"Output policy {name!r} must be an OutputPolicy or mapping")
        normalized[name] = policy
    return MappingProxyType(normalized)


def _deep_freeze(value: Any) -> Any:
    """Copy nested containers into immutable equivalents."""

    if isinstance(value, Mapping):
        return MappingProxyType({key: _deep_freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_deep_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_deep_freeze(item) for item in value)
    return value


def _choice(*choices: str) -> ParameterValidator:
    def validate(value: Any) -> str:
        if value not in choices:
            raise ValueError(f"must be one of: {', '.join(choices)}")
        return value

    return validate


def _optional_choice(*choices: str) -> ParameterValidator:
    """Accept one of ``choices`` or ``None`` for an opt-in upstream setting."""

    def validate(value: Any) -> str | None:
        if value is None:
            return None
        if value not in choices:
            raise ValueError(f"must be null or one of: {', '.join(choices)}")
        return value

    return validate


def _annotation_tools(value: Any) -> tuple[str, ...]:
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        raise ValueError("must be a list of annotation tool names")
    tools = tuple(str(item) for item in value)
    unknown = sorted(set(tools) - set(ANNOTATION_TOOLS))
    if unknown:
        raise ValueError(f"must contain only: {', '.join(ANNOTATION_TOOLS)}")
    if len(set(tools)) != len(tools):
        raise ValueError("must not repeat a tool")
    return tools


def _optional_positive_int(value: Any) -> int | None:
    if value is None:
        return None
    return _positive_int(value)


def _optional_string(value: Any) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError("must be a non-empty string or null")
    return value


def _optional_path(value: Any) -> Path | None:
    if value is None:
        return None
    if not isinstance(value, (str, Path)):
        raise ValueError("must be a filesystem path or null")
    return Path(value)


def _identifier(value: Any) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", value):
        raise ValueError("must be a simple column name")
    return value


def _optional_identifier(value: Any) -> str | None:
    if value is None:
        return None
    return _identifier(value)


def _factor_level(value: Any) -> str:
    if not isinstance(value, str) or not re.fullmatch(
        r"[A-Za-z0-9][A-Za-z0-9_.-]*", value
    ):
        raise ValueError("must be a simple factor level")
    return value


def _positive_int(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("must be a positive integer")
    return value


def _non_negative_int(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("must be a non-negative integer")
    return value


def _positive_number(value: Any) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("must be a positive number")
    if not math.isfinite(value) or value <= 0:
        raise ValueError("must be a positive finite number")
    return value


def _non_negative_number(value: Any) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("must be a non-negative number")
    if not math.isfinite(value) or value < 0:
        raise ValueError("must be a non-negative finite number")
    return value


def _percentage(value: Any) -> int | float:
    value = _non_negative_number(value)
    if value > 100:
        raise ValueError("must be between 0 and 100")
    return value


def _probability(value: Any) -> int | float:
    value = _positive_number(value)
    if value > 1:
        raise ValueError("must be in (0, 1]")
    return value


def _boolean(value: Any) -> bool:
    if not isinstance(value, bool):
        raise ValueError("must be a boolean")
    return value


def _covariates(value: Any) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError("must be a list of simple column names")
    return tuple(_identifier(item) for item in value)


_PATH_OPTIONS = {
    "params_file": _optional_path,
    "nextflow_config": _optional_path,
}
_IMAGE_OPTION = {"image": _optional_string}


STAGE_SPECS: tuple[StageSpec, ...] = (
    StageSpec(
        id="bulk.rnaseq",
        modality="bulk",
        maturity="selected",
        inputs={"samplesheet": "bulk.rnaseq-samplesheet"},
        outputs={
            "results": "bulk.rnaseq-results",
            "counts": "bulk.gene-counts",
            "multiqc_report": "qc.multiqc-report",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "counts": OutputPolicy("file", non_empty=True),
            "multiqc_report": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "params_file": None,
            "nextflow_config": None,
            "pseudo_aligner": None,
            "skip_alignment": False,
            "salmon_index": None,
        },
        validators={
            **_PATH_OPTIONS,
            "pseudo_aligner": _optional_choice(*PSEUDO_ALIGNERS),
            "skip_alignment": _boolean,
            "salmon_index": _optional_path,
        },
        required_executables=("nextflow",),
        command_factory=bulk_rnaseq_command,
        supports_resume=True,
    ),
    StageSpec(
        id="bulk.salmon",
        modality="bulk",
        maturity="ready",
        # A weaker contract than bulk.rnaseq-samplesheet: selective alignment
        # infers strandedness, so no strandedness column is required. Typing it
        # separately stops a salmon sheet being wired into bulk.rnaseq, where
        # nf-core would reject it.
        inputs={"samplesheet": "bulk.fastq-samplesheet"},
        outputs={
            "results": "bulk.salmon-results",
            "counts": "bulk.gene-counts",
            "tx_counts": "bulk.transcript-counts",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "counts": OutputPolicy("file", non_empty=True),
            "tx_counts": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "fasta": None,
            "gtf": None,
            "salmon_index": None,
            "tx2gene": None,
            "image": None,
            "libtype": "A",
            "kmer_len": 31,
            "gencode": False,
            "nextflow_config": None,
        },
        validators={
            "fasta": _optional_path,
            "gtf": _optional_path,
            "salmon_index": _optional_path,
            "tx2gene": _optional_path,
            "image": _optional_string,
            "libtype": _choice(*SALMON_LIBTYPES),
            "kmer_len": _positive_int,
            "gencode": _boolean,
            "nextflow_config": _optional_path,
        },
        required_executables=("nextflow",),
        required_images=("images.salmon",),
        command_factory=bulk_salmon_command,
        supports_resume=True,
    ),
    StageSpec(
        id="bulk.align",
        modality="bulk",
        maturity="ready",
        inputs={
            "samplesheet": "bulk.fastq-samplesheet",
            "star_index": "bulk.star-index",
        },
        outputs={
            "results": "bulk.align-results",
            "alignments": "bulk.alignment-bams",
            "counts": "bulk.gene-counts",
            "strandedness": "bulk.strandedness-report",
            "logs": "bulk.alignment-logs",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "alignments": OutputPolicy("directory", non_empty=True),
            "counts": OutputPolicy("file", non_empty=True),
            "strandedness": OutputPolicy("file", non_empty=True),
            "logs": OutputPolicy("directory", non_empty=True),
        },
        defaults={
            "image": None,
            "strandedness": "auto",
            "two_pass": True,
            "threads": 4,
            "memory_gb": 32,
            "nextflow_config": None,
        },
        validators={
            "image": _optional_string,
            "strandedness": _choice(*STRANDEDNESS),
            "two_pass": _boolean,
            "threads": _positive_int,
            "memory_gb": _positive_int,
            "nextflow_config": _optional_path,
        },
        required_executables=("nextflow",),
        required_images=("images.star",),
        command_factory=bulk_align_command,
        supports_resume=True,
    ),
    StageSpec(
        id="bulk.rnavar",
        modality="bulk",
        maturity="selected",
        inputs={"samplesheet": "bulk.rnavar-samplesheet"},
        # Supplying a prebuilt reference is optional: rnavar derives whatever is
        # missing, so a project can wire ${ref.star_index} from
        # bulk.star-reference to avoid rebuilding it on every run, or omit it.
        optional_inputs={
            "star_index": "bulk.star-index",
            "fasta_fai": "reference.fasta-index",
            "dict": "reference.sequence-dictionary",
        },
        outputs={
            "results": "bulk.rnavar-results",
            "variants": "bulk.variant-calls",
            "multiqc_report": "qc.multiqc-report",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "variants": OutputPolicy("directory", non_empty=True),
            "multiqc_report": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "genome": None,
            "fasta": None,
            "gtf": None,
            "dbsnp": None,
            "known_indels": None,
            "skip_baserecalibration": False,
            "tools": (),
            "snpeff_cache": None,
            "vep_cache": None,
            "generate_gvcf": False,
            "params_file": None,
            "nextflow_config": None,
        },
        validators={
            "genome": _optional_string,
            "fasta": _optional_path,
            "gtf": _optional_path,
            "dbsnp": _optional_path,
            "known_indels": _optional_path,
            "skip_baserecalibration": _boolean,
            "tools": _annotation_tools,
            "snpeff_cache": _optional_path,
            "vep_cache": _optional_path,
            "generate_gvcf": _boolean,
            "params_file": _optional_path,
            "nextflow_config": _optional_path,
        },
        required_executables=("nextflow",),
        command_factory=bulk_rnavar_command,
        supports_resume=True,
    ),
    StageSpec(
        id="bulk.star-reference",
        modality="bulk",
        maturity="ready",
        inputs={
            "fasta": "reference.genome-fasta",
            "gtf": "reference.annotation-gtf",
        },
        outputs={
            "results": "bulk.star-reference-results",
            "star_index": "bulk.star-index",
            "fasta_fai": "reference.fasta-index",
            "dict": "reference.sequence-dictionary",
            "manifest": "reference.index-manifest",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "star_index": OutputPolicy("directory", non_empty=True),
            "fasta_fai": OutputPolicy("file", non_empty=True),
            "dict": OutputPolicy("file", non_empty=True),
            "manifest": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "image": None,
            "read_length": 100,
            "sjdb_overhang": None,
            "threads": 4,
            "memory_gb": 32,
            "genome_sa_index_nbases": None,
            "nextflow_config": None,
        },
        validators={
            "image": _optional_string,
            "read_length": _positive_int,
            "sjdb_overhang": _optional_positive_int,
            "threads": _positive_int,
            "memory_gb": _positive_int,
            "genome_sa_index_nbases": _optional_positive_int,
            "nextflow_config": _optional_path,
        },
        required_executables=("nextflow",),
        required_images=("images.star",),
        command_factory=bulk_star_reference_command,
        supports_resume=True,
    ),
    StageSpec(
        id="bulk.de",
        modality="bulk",
        maturity="ready",
        inputs={"counts": "bulk.gene-counts", "metadata": "sample.metadata"},
        outputs={
            "results": "bulk.differential-expression-results",
            "de_results": "bulk.differential-expression-table",
            "contrast_index": "bulk.differential-expression-index",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "de_results": OutputPolicy("file", non_empty=True),
            "contrast_index": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "method": "deseq2",
            "image": None,
            "covariates": (),
            "padj": 0.05,
            "lfc": 1.0,
            "top_genes": 50,
            "contrasts": "single",
        },
        validators={
            "method": _choice("deseq2", "edger", "limma"),
            "image": _optional_string,
            "design": _identifier,
            "reference": _factor_level,
            "test": _factor_level,
            "covariates": _covariates,
            "padj": _probability,
            "lfc": _non_negative_number,
            "top_genes": _positive_int,
            "contrasts": _choice(*CONTRAST_MODES),
        },
        required_executables=("docker",),
        required_images=("images.bulk_r",),
        command_factory=bulk_de_command,
    ),
    StageSpec(
        id="bulk.enrichment",
        modality="bulk",
        maturity="ready",
        inputs={
            "de_results": "bulk.differential-expression-table",
            "genesets": "gene-sets.gmt",
        },
        outputs={"results": "bulk.gene-set-enrichment"},
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
        },
        defaults={
            "image": None,
            "mode": "ora",
            "padj": 0.05,
            "lfc": 1.0,
            "min_size": 10,
            "max_size": 500,
            "adjust": "BH",
        },
        validators={
            "image": _optional_string,
            "mode": _choice("ora", "gsea"),
            "padj": _probability,
            "lfc": _non_negative_number,
            "min_size": _positive_int,
            "max_size": _positive_int,
            "adjust": _choice(
                "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
            ),
        },
        required_executables=("docker",),
        required_images=("images.bulk_r",),
        command_factory=bulk_enrichment_command,
    ),
    StageSpec(
        id="single-cell.scrnaseq",
        modality="single-cell",
        maturity="selected",
        inputs={"samplesheet": "single-cell.scrnaseq-samplesheet"},
        outputs={
            "results": "single-cell.scrnaseq-results",
            "matrix": "single-cell.matrix",
            "multiqc_report": "qc.multiqc-report",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "matrix": OutputPolicy("file", non_empty=True),
            "multiqc_report": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "aligner": "simpleaf",
            "protocol": None,
            "params_file": None,
            "nextflow_config": None,
            "simpleaf_index": None,
            "txp2gene": None,
        },
        validators={
            "aligner": _choice("simpleaf", "star", "cellranger"),
            "protocol": _optional_string,
            "simpleaf_index": _optional_path,
            "txp2gene": _optional_path,
            **_PATH_OPTIONS,
        },
        required_executables=("nextflow",),
        command_factory=scrnaseq_command,
        supports_resume=True,
    ),
    StageSpec(
        id="single-cell.alevin",
        modality="single-cell",
        maturity="ready",
        inputs={"samplesheet": "single-cell.scrnaseq-samplesheet"},
        outputs={
            "results": "single-cell.alevin-results",
            "matrix": "single-cell.matrix",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "matrix": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "fasta": None,
            "gtf": None,
            "simpleaf_index": None,
            "whitelist": None,
            "image": None,
            "single_cell_image": None,
            "chemistry": "10xv3",
            "resolution": "cr-like",
            "rlen": 91,
            "nextflow_config": None,
        },
        validators={
            "fasta": _optional_path,
            "gtf": _optional_path,
            "simpleaf_index": _optional_path,
            "whitelist": _optional_path,
            "image": _optional_string,
            "single_cell_image": _optional_string,
            "chemistry": _choice(*ALEVIN_CHEMISTRIES),
            "resolution": _choice(*ALEVIN_RESOLUTIONS),
            "rlen": _positive_int,
            "nextflow_config": _optional_path,
        },
        required_executables=("nextflow",),
        required_images=("images.salmon", "images.single_cell_python"),
        command_factory=alevin_command,
        supports_resume=True,
    ),
    StageSpec(
        id="single-cell.scanpy",
        modality="single-cell",
        maturity="ready",
        inputs={"input": "single-cell.matrix"},
        # Cell metadata is optional because the batch column may already be
        # present in the matrix's obs; analysis_command mounts it only when given.
        optional_inputs={"metadata": "single-cell.cell-metadata"},
        outputs={"h5ad": "single-cell.h5ad"},
        output_policies={"h5ad": OutputPolicy("file", non_empty=True)},
        defaults={
            "image": None,
            "min_genes": 200,
            "min_cells": 3,
            "max_mito_pct": 20.0,
            "resolution": 1.0,
            "integration": "none",
            "batch_column": None,
            "barcode_column": "barcode",
        },
        validators={
            **_IMAGE_OPTION,
            "min_genes": _non_negative_int,
            "min_cells": _non_negative_int,
            "max_mito_pct": _percentage,
            "resolution": _positive_number,
            "integration": _choice(*INTEGRATION_METHODS),
            "batch_column": _optional_identifier,
            "barcode_column": _identifier,
        },
        required_executables=("docker",),
        required_images=("images.single_cell_python",),
        command_factory=scanpy_command,
    ),
    StageSpec(
        id="single-cell.pseudobulk",
        modality="single-cell",
        maturity="ready",
        inputs={"h5ad": "single-cell.h5ad"},
        outputs={
            "counts": "bulk.gene-counts",
            "metadata": "sample.metadata",
        },
        output_policies={
            "counts": OutputPolicy("file", non_empty=True),
            "metadata": OutputPolicy("file", non_empty=True),
        },
        defaults={"image": None},
        validators={
            **_IMAGE_OPTION,
            "sample_column": _identifier,
            "design": _identifier,
        },
        required_executables=("docker",),
        required_images=("images.single_cell_python",),
        command_factory=pseudobulk_command,
    ),
    StageSpec(
        id="single-cell.pseudobulk-de",
        modality="single-cell",
        maturity="ready",
        inputs={"h5ad": "single-cell.h5ad"},
        outputs={
            "results": "single-cell.pseudobulk-de-results",
            "de_results": "bulk.differential-expression-table",
        },
        output_policies={
            "results": OutputPolicy("directory", non_empty=True),
            "de_results": OutputPolicy("file", non_empty=True),
        },
        defaults={
            "single_cell_image": None,
            "bulk_image": None,
            "padj": 0.05,
            "lfc": 1.0,
            "top_genes": 50,
            "nextflow_config": None,
        },
        validators={
            "sample_column": _identifier,
            "design": _identifier,
            "reference": _factor_level,
            "test": _factor_level,
            "single_cell_image": _optional_string,
            "bulk_image": _optional_string,
            "padj": _probability,
            "lfc": _non_negative_number,
            "top_genes": _positive_int,
            "nextflow_config": _optional_path,
        },
        required_executables=("nextflow",),
        required_images=("images.single_cell_python", "images.bulk_r"),
        command_factory=pseudobulk_de_command,
        supports_resume=True,
    ),
)

_STAGE_REGISTRY = MappingProxyType({spec.uses: spec for spec in STAGE_SPECS})


def get_stage_spec(uses: str) -> StageSpec:
    """Return the exact registered stage, failing on unknown ``uses`` keys."""

    try:
        return _STAGE_REGISTRY[uses]
    except KeyError as exc:
        known = ", ".join(_STAGE_REGISTRY)
        raise TxSuiteError(f"Unknown project stage {uses!r}; known stages: {known}") from exc


def list_stage_specs(modality: str | None = None) -> tuple[StageSpec, ...]:
    """List stages in stable registry order, optionally filtered by modality."""

    if modality is None:
        return STAGE_SPECS
    return tuple(spec for spec in STAGE_SPECS if spec.modality == modality)


def validate_stage_parameters(
    spec: StageSpec | str, params: Mapping[str, Any] | None
) -> dict[str, Any]:
    """Apply defaults and validators, rejecting unknown or missing parameters."""

    if isinstance(spec, str):
        spec = get_stage_spec(spec)
    if params is None:
        params = {}
    if not isinstance(params, Mapping):
        raise TxSuiteError(f"Parameters for stage {spec.uses} must be a mapping")
    allowed = set(spec.defaults) | set(spec.validators)
    unknown = sorted(set(params) - allowed)
    if unknown:
        raise TxSuiteError(
            f"Unknown parameter(s) for stage {spec.uses}: {', '.join(unknown)}"
        )
    required = set(spec.validators) - set(spec.defaults)
    missing = sorted(required - set(params))
    if missing:
        raise TxSuiteError(
            f"Missing required parameter(s) for stage {spec.uses}: {', '.join(missing)}"
        )
    normalized = dict(spec.defaults)
    normalized.update(params)
    for name, validator in spec.validators.items():
        if name not in normalized:
            continue
        try:
            normalized[name] = validator(normalized[name])
        except (TypeError, ValueError) as exc:
            raise TxSuiteError(
                f"Invalid parameter {name!r} for stage {spec.uses}: {exc}"
            ) from exc
    if "reference" in normalized and "test" in normalized:
        if normalized["reference"] == normalized["test"]:
            raise TxSuiteError(
                f"Invalid parameters for stage {spec.uses}: reference and test must differ"
            )
    if "min_size" in normalized and "max_size" in normalized:
        if normalized["max_size"] < normalized["min_size"]:
            raise TxSuiteError(
                f"Invalid parameters for stage {spec.uses}: max_size must be >= min_size"
            )
    return normalized


__all__ = [
    "OutputPolicy",
    "ParameterValidator",
    "STAGE_SPECS",
    "StageCommandFactory",
    "StageSpec",
    "get_stage_spec",
    "list_stage_specs",
    "validate_stage_parameters",
]
