"""RNA-seq short variant discovery through the pinned nf-core/rnavar release.

TxSuite launches the upstream pipeline rather than reimplementing the GATK
best-practice chain (STAR two-pass, MarkDuplicates, SplitNCigarReads, BQSR,
HaplotypeCaller, VariantFiltration). What this module owns is the part that
fails cheaply: rejecting reference, recalibration, and annotation combinations
that would otherwise die hours into a run.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError

ANNOTATION_TOOLS = ("snpeff", "vep")


def validate_rnavar_samplesheet(path: Path) -> int:
    """Validate the FASTQ-entry samplesheet the rnavar stage accepts.

    rnavar also supports BAM, CRAM, and VCF entry points. TxSuite models only
    the FASTQ one for now, so a samplesheet built for another entry point is
    rejected here instead of failing later inside the pipeline.
    """

    if not path.is_file():
        raise TxSuiteError(f"Samplesheet does not exist: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or ())
        missing = {"sample", "fastq_1"} - fieldnames
        if missing:
            raise TxSuiteError(
                f"Samplesheet is missing columns: {', '.join(sorted(missing))}"
            )
        unsupported = sorted(fieldnames & {"bam", "cram", "vcf"})
        if unsupported:
            raise TxSuiteError(
                "TxSuite supports only the FASTQ entry point; samplesheet declares: "
                + ", ".join(unsupported)
            )
        seen: set[str] = set()
        rows = 0
        for line, row in enumerate(reader, start=2):
            sample = (row.get("sample") or "").strip()
            if not sample or not (row.get("fastq_1") or "").strip():
                raise TxSuiteError(
                    f"Samplesheet line {line} requires sample and fastq_1"
                )
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sample):
                raise TxSuiteError(
                    f"Samplesheet line {line} sample may contain only letters, "
                    "numbers, ., _ and -"
                )
            if sample in seen:
                raise TxSuiteError(f"Samplesheet line {line} repeats sample {sample!r}")
            seen.add(sample)
            rows += 1
    if not rows:
        raise TxSuiteError("Samplesheet has no data rows")
    return rows


def rnavar_workflow_command(
    config: dict[str, Any],
    *,
    samplesheet: Path,
    outdir: Path,
    genome: str | None = None,
    fasta: Path | None = None,
    gtf: Path | None = None,
    star_index: Path | None = None,
    fasta_fai: Path | None = None,
    sequence_dictionary: Path | None = None,
    dbsnp: Path | None = None,
    known_indels: Path | None = None,
    skip_baserecalibration: bool = False,
    tools: tuple[str, ...] = (),
    snpeff_cache: Path | None = None,
    vep_cache: Path | None = None,
    generate_gvcf: bool = False,
    params_file: Path | None = None,
    nextflow_config: Path | None = None,
    resume: bool = False,
    check_inputs: bool = True,
) -> list[str]:
    """Build the pinned nf-core/rnavar command.

    Reference selection is exclusive: either an iGenomes ``genome`` key, or a
    ``fasta`` and ``gtf`` pair. Base quality score recalibration needs known
    sites, and each annotation tool needs its cache; both are rejected here so
    the failure arrives before alignment rather than after it.
    """

    if bool(genome) == bool(fasta or gtf):
        raise TxSuiteError(
            "Pass either --genome, or both --fasta and --gtf, to select a reference"
        )
    if bool(fasta) != bool(gtf):
        raise TxSuiteError("A custom reference requires both --fasta and --gtf")
    if genome is not None and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", genome):
        raise TxSuiteError("Genome key may contain only letters, numbers, ., _ and -")

    unknown_tools = sorted(set(tools) - set(ANNOTATION_TOOLS))
    if unknown_tools:
        raise TxSuiteError(
            f"Unknown annotation tool(s): {', '.join(unknown_tools)}; "
            f"supported: {', '.join(ANNOTATION_TOOLS)}"
        )
    if len(set(tools)) != len(tools):
        raise TxSuiteError("Annotation tools must be unique")
    if "snpeff" in tools and snpeff_cache is None:
        raise TxSuiteError("snpEff annotation requires --snpeff-cache")
    if "vep" in tools and vep_cache is None:
        raise TxSuiteError("VEP annotation requires --vep-cache")

    if not skip_baserecalibration and dbsnp is None and known_indels is None:
        raise TxSuiteError(
            "Base recalibration requires --dbsnp or --known-indels; pass "
            "--skip-baserecalibration to run without it"
        )
    if skip_baserecalibration and known_indels is not None:
        # dbSNP is deliberately still allowed here: rnavar feeds it to
        # HaplotypeCaller for rsID annotation, which is independent of BQSR.
        # Known indels have no consumer once recalibration is skipped.
        raise TxSuiteError(
            "Known indels are unused when base recalibration is skipped"
        )

    if check_inputs:
        validate_rnavar_samplesheet(samplesheet)
        for label, path, is_dir in (
            ("Genome FASTA", fasta, False),
            ("Annotation GTF", gtf, False),
            ("STAR index", star_index, True),
            ("FASTA index", fasta_fai, False),
            ("Sequence dictionary", sequence_dictionary, False),
            ("dbSNP VCF", dbsnp, False),
            ("Known indels VCF", known_indels, False),
            ("snpEff cache", snpeff_cache, True),
            ("VEP cache", vep_cache, True),
        ):
            if path is None:
                continue
            if is_dir and not path.is_dir():
                raise TxSuiteError(f"{label} does not exist: {path}")
            if not is_dir and not path.is_file():
                raise TxSuiteError(f"{label} does not exist: {path}")
    if params_file is not None and not params_file.is_file():
        raise TxSuiteError(f"Params file does not exist: {params_file}")
    if nextflow_config is not None and not nextflow_config.is_file():
        raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")

    pipeline = config["pipelines"]["variants"]
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
    if genome is not None:
        command.extend(["--genome", genome])
    else:
        command.extend(
            ["--fasta", str(fasta.resolve()), "--gtf", str(gtf.resolve())]
        )
    for flag, path in (
        ("--star_index", star_index),
        ("--fasta_fai", fasta_fai),
        ("--dict", sequence_dictionary),
        ("--dbsnp", dbsnp),
        ("--known_indels", known_indels),
        ("--snpeff_cache", snpeff_cache),
        ("--vep_cache", vep_cache),
    ):
        if path is not None:
            command.extend([flag, str(path.resolve())])
    if skip_baserecalibration:
        command.append("--skip_baserecalibration")
    if tools:
        command.extend(["--tools", ",".join(tools)])
    if generate_gvcf:
        command.append("--generate_gvcf")
    if params_file is not None:
        command.extend(["-params-file", str(params_file.resolve())])
    if nextflow_config is not None:
        command.extend(["-c", str(nextflow_config.resolve())])
    if resume:
        command.append("-resume")
    return command


__all__ = [
    "ANNOTATION_TOOLS",
    "rnavar_workflow_command",
    "validate_rnavar_samplesheet",
]
