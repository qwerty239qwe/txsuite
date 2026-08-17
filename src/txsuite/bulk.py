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
PSEUDO_ALIGNERS = ("salmon", "kallisto")
CONTRAST_MODES = ("single", "vs-reference", "all-pairs")


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
    pseudo_aligner: str | None = None,
    skip_alignment: bool = False,
    salmon_index: Path | None = None,
    resume: bool = False,
) -> list[str]:
    if pseudo_aligner is not None and pseudo_aligner not in PSEUDO_ALIGNERS:
        raise TxSuiteError(
            f"Pseudo-aligner must be one of: {', '.join(PSEUDO_ALIGNERS)}"
        )
    if skip_alignment and pseudo_aligner is None:
        raise TxSuiteError("Skipping alignment requires a pseudo-aligner")
    if salmon_index is not None:
        if pseudo_aligner != "salmon":
            raise TxSuiteError("A salmon index requires --pseudo-aligner salmon")
        if not salmon_index.is_dir():
            raise TxSuiteError(f"Salmon index does not exist: {salmon_index}")
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
    if pseudo_aligner is not None:
        command.extend(["--pseudo_aligner", pseudo_aligner])
    if skip_alignment:
        command.append("--skip_alignment")
    if salmon_index is not None:
        command.extend(["--salmon_index", str(salmon_index.resolve())])
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


SALMON_LIBTYPES = (
    "A",
    "IU",
    "ISF",
    "ISR",
    "MU",
    "MSF",
    "MSR",
    "OU",
    "OSF",
    "OSR",
    "U",
    "SF",
    "SR",
)


def validate_quant_samplesheet(path: Path) -> int:
    """Validate the reduced samplesheet the native quantification DAGs accept.

    Selective alignment infers strandedness from ``--libType A``, so unlike
    :func:`validate_samplesheet` this contract requires only ``sample`` and
    ``fastq_1``. Reading the file here means a malformed samplesheet fails
    before Nextflow starts rather than inside a process.
    """

    if not path.is_file():
        raise TxSuiteError(f"Samplesheet does not exist: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"sample", "fastq_1"} - set(reader.fieldnames or ())
        if missing:
            raise TxSuiteError(
                f"Samplesheet is missing columns: {', '.join(sorted(missing))}"
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


def salmon_workflow_command(
    config: dict[str, Any],
    *,
    samplesheet: Path,
    outdir: Path,
    fasta: Path | None = None,
    gtf: Path | None = None,
    salmon_index: Path | None = None,
    tx2gene: Path | None = None,
    salmon_image: str | None = None,
    libtype: str = "A",
    kmer_len: int = 31,
    gencode: bool = False,
    nextflow_config: Path | None = None,
    resume: bool = False,
    check_inputs: bool = True,
) -> list[str]:
    """Build the native salmon selective-alignment Nextflow DAG.

    Either ``salmon_index`` or both ``fasta`` and ``gtf`` must be supplied. A
    prebuilt index carries no annotation, so it also requires ``tx2gene``;
    when the index is built here the table is derived from the GTF.
    """

    if libtype not in SALMON_LIBTYPES:
        raise TxSuiteError(f"Salmon library type must be one of: {', '.join(SALMON_LIBTYPES)}")
    if kmer_len < 1 or kmer_len % 2 == 0:
        raise TxSuiteError("Salmon k-mer length must be a positive odd number")
    if bool(salmon_index) == bool(fasta or gtf):
        raise TxSuiteError(
            "Pass either --salmon-index, or both --fasta and --gtf to build one"
        )
    if bool(fasta) != bool(gtf):
        raise TxSuiteError("Building a salmon index requires both --fasta and --gtf")
    if salmon_index is not None and tx2gene is None:
        raise TxSuiteError(
            "A prebuilt salmon index also requires --tx2gene for gene-level counts"
        )
    if check_inputs:
        validate_quant_samplesheet(samplesheet)
        for label, path, is_dir in (
            ("Salmon index", salmon_index, True),
            ("Genome FASTA", fasta, False),
            ("Annotation GTF", gtf, False),
            ("Transcript-to-gene table", tx2gene, False),
        ):
            if path is None:
                continue
            if is_dir and not path.is_dir():
                raise TxSuiteError(f"{label} does not exist: {path}")
            if not is_dir and not path.is_file():
                raise TxSuiteError(f"{label} does not exist: {path}")
    if nextflow_config is not None and not nextflow_config.is_file():
        raise TxSuiteError(f"Nextflow config does not exist: {nextflow_config}")

    image = salmon_image or config["images"]["salmon"]
    if not image.strip():
        raise TxSuiteError("Workflow images cannot be empty")

    workflow = resources.files("txsuite.resources.nextflow").joinpath("bulk_salmon.nf")
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
            "--samplesheet",
            str(samplesheet.resolve()),
            "--outdir",
            str(outdir.resolve()),
            "--salmon_image",
            image,
            "--salmon_libtype",
            libtype,
            "--salmon_kmer_len",
            str(kmer_len),
            "--salmon_gencode",
            "true" if gencode else "false",
        ]
    )
    if salmon_index is not None:
        command.extend(["--salmon_index", str(salmon_index.resolve())])
    else:
        command.extend(
            ["--fasta", str(fasta.resolve()), "--gtf", str(gtf.resolve())]
        )
    if tx2gene is not None:
        command.extend(["--tx2gene", str(tx2gene.resolve())])
    return command


def deseq2_command(
    *,
    image: str,
    counts: Path,
    metadata: Path,
    outdir: Path,
    design: str | None = None,
    reference: str | None = None,
    test: str | None = None,
    formula: str | None = None,
    coefficient: str | None = None,
    covariates: tuple[str, ...] = (),
    padj: float = 0.05,
    lfc: float = 1.0,
    top_genes: int = 50,
    contrasts: str = "single",
    check_inputs: bool = True,
) -> list[str]:
    if check_inputs:
        for label, path in (("Counts", counts), ("Metadata", metadata)):
            if not path.is_file():
                raise TxSuiteError(f"{label} file does not exist: {path}")
    if contrasts not in CONTRAST_MODES:
        raise TxSuiteError(f"Contrasts must be one of: {', '.join(CONTRAST_MODES)}")
    advanced = formula is not None or coefficient is not None
    if advanced and contrasts != "single":
        raise TxSuiteError(
            "Formula mode has no design levels to expand; contrasts must be 'single'"
        )
    if advanced:
        if not formula or not coefficient:
            raise TxSuiteError("Formula and coefficient must be used together")
        if any(value is not None for value in (design, reference, test)) or covariates:
            raise TxSuiteError(
                "Formula mode cannot be combined with design, contrast levels, or covariates"
            )
        if not re.fullmatch(r"~[A-Za-z0-9_.+*: ]+", formula) or not re.search(
            r"[A-Za-z]", formula
        ):
            raise TxSuiteError("Formula supports column names and +, *, :, 0, or 1")
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.:]*", coefficient):
            raise TxSuiteError("Coefficient must be a simple model coefficient name")
    else:
        if design is None or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", design):
            raise TxSuiteError("Design must be a simple metadata column name")
        if not reference or not test or reference == test:
            raise TxSuiteError(
                "Reference and test levels must be non-empty and different"
            )
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
        design or "",
        reference or "",
        test or "",
        "/output",
        str(padj),
        str(lfc),
        str(top_genes),
        ",".join(covariates),
        formula or "",
        coefficient or "",
        contrasts,
    ]


def differential_expression_command(
    *,
    method: str,
    image: str,
    counts: Path,
    metadata: Path,
    outdir: Path,
    design: str | None = None,
    reference: str | None = None,
    test: str | None = None,
    formula: str | None = None,
    coefficient: str | None = None,
    covariates: tuple[str, ...] = (),
    padj: float = 0.05,
    lfc: float = 1.0,
    top_genes: int = 50,
    contrasts: str = "single",
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
        formula=formula,
        coefficient=coefficient,
        outdir=outdir,
        covariates=covariates,
        padj=padj,
        lfc=lfc,
        top_genes=top_genes,
        contrasts=contrasts,
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
    check_inputs: bool = True,
) -> list[str]:
    if check_inputs:
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


def build_salmon_image(tag: str, *, run_dir: Path) -> None:
    """Build the quantification image shared by the bulk and single-cell DAGs."""

    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    package = resources.files("txsuite.resources.salmon")
    with tempfile.TemporaryDirectory(prefix="txsuite-salmon-") as directory:
        context = Path(directory)
        for name in ("Dockerfile", "merge_quants.py"):
            (context / name).write_text(
                package.joinpath(name).read_text(encoding="utf-8"), encoding="utf-8"
            )
        run_command(
            ["docker", "build", "--tag", tag, str(context)],
            run_dir=run_dir,
            task="env.build.salmon",
            backend="docker",
            inputs={
                "base": (
                    "condaforge/miniforge3:26.3.2-3@"
                    "sha256:532f6ee7a858b009dc895f8313eb6ed875f05a455ec57d832aa6b4c66e2799b9"
                ),
                "salmon": "2.5.1",
                "simpleaf": "0.28.0",
                "alevin-fry": "0.18.0",
                "piscem": "0.22.1",
                "gffread": "0.12.9",
            },
            outputs={"image": tag},
            artifacts=[{"kind": "container-image", "label": "salmon", "path": tag}],
        )


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
