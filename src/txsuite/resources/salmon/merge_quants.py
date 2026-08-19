"""Merge per-sample salmon quantifications into transcript and gene matrices.

Salmon writes one ``quant.sf`` per sample. Downstream TxSuite stages consume a
single matrix whose columns are samples, so this script performs the join and
the transcript-to-gene summarization in one dependency-free pass. Counts are
summed per gene and rounded to integers because ``bulk.de`` backends (DESeq2,
edgeR) require an integer count matrix; TPM is summed without rounding.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


QUANT_FILE = "quant.sf"


class MergeError(RuntimeError):
    """A quantification input is missing, malformed, or mutually inconsistent."""


def read_tx2gene(path: Path) -> dict[str, str]:
    """Read a two-column transcript-to-gene table, tolerating a header row."""

    mapping: dict[str, str] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise MergeError(
                    f"{path}:{line_number}: expected transcript and gene columns"
                )
            transcript, gene = row[0].strip(), row[1].strip()
            if not transcript or not gene:
                continue
            if transcript.lower() in {"transcript_id", "tx", "txname", "name"}:
                continue
            existing = mapping.get(transcript)
            if existing is not None and existing != gene:
                raise MergeError(
                    f"{path}:{line_number}: transcript {transcript!r} maps to both "
                    f"{existing!r} and {gene!r}"
                )
            mapping[transcript] = gene
    if not mapping:
        raise MergeError(f"Transcript-to-gene table is empty: {path}")
    return mapping


def read_quant(path: Path) -> dict[str, tuple[float, float]]:
    """Read one ``quant.sf`` into ``{transcript: (num_reads, tpm)}``."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"Name", "TPM", "NumReads"} - set(reader.fieldnames or ())
        if missing:
            raise MergeError(
                f"{path}: missing salmon columns: {', '.join(sorted(missing))}"
            )
        values: dict[str, tuple[float, float]] = {}
        for line_number, row in enumerate(reader, start=2):
            name = (row["Name"] or "").strip()
            if not name:
                continue
            try:
                values[name] = (float(row["NumReads"]), float(row["TPM"]))
            except (TypeError, ValueError) as exc:
                raise MergeError(f"{path}:{line_number}: non-numeric salmon value") from exc
    if not values:
        raise MergeError(f"No quantified transcripts in {path}")
    return values


def _quant_path(quant_dir: Path) -> Path:
    candidate = quant_dir / QUANT_FILE
    if candidate.is_file():
        return candidate
    if quant_dir.is_file() and quant_dir.name == QUANT_FILE:
        return quant_dir
    raise MergeError(f"No {QUANT_FILE} under {quant_dir}")


def _write_matrix(
    path: Path,
    row_label: str,
    samples: list[str],
    rows: dict[str, dict[str, float]],
    *,
    integer: bool,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_label, *samples])
        for key in sorted(rows):
            per_sample = rows[key]
            if integer:
                cells = [str(int(round(per_sample.get(sample, 0.0)))) for sample in samples]
            else:
                cells = [f"{per_sample.get(sample, 0.0):.6f}" for sample in samples]
            writer.writerow([key, *cells])


def merge(quant_dirs: list[Path], tx2gene_path: Path, outdir: Path) -> dict[str, Path]:
    """Write transcript and gene matrices, returning the paths by artifact name."""

    tx2gene = read_tx2gene(tx2gene_path)
    samples: list[str] = []
    transcript_counts: dict[str, dict[str, float]] = {}
    transcript_tpm: dict[str, dict[str, float]] = {}
    gene_counts: dict[str, dict[str, float]] = {}
    gene_tpm: dict[str, dict[str, float]] = {}
    unmapped: set[str] = set()

    for quant_dir in quant_dirs:
        sample = quant_dir.name
        if sample in samples:
            raise MergeError(f"Duplicate sample name: {sample}")
        samples.append(sample)
        for transcript, (reads, tpm) in read_quant(_quant_path(quant_dir)).items():
            transcript_counts.setdefault(transcript, {})[sample] = reads
            transcript_tpm.setdefault(transcript, {})[sample] = tpm
            gene = tx2gene.get(transcript)
            if gene is None:
                unmapped.add(transcript)
                continue
            gene_counts.setdefault(gene, {})[sample] = (
                gene_counts.get(gene, {}).get(sample, 0.0) + reads
            )
            gene_tpm.setdefault(gene, {})[sample] = (
                gene_tpm.get(gene, {}).get(sample, 0.0) + tpm
            )

    if not gene_counts:
        raise MergeError(
            "No quantified transcript matched the transcript-to-gene table; the "
            "index and annotation are probably from different sources"
        )
    if unmapped:
        print(
            f"warning: {len(unmapped)} transcripts had no gene assignment and were "
            "excluded from the gene matrices",
            file=sys.stderr,
        )

    outdir.mkdir(parents=True, exist_ok=True)
    written = {
        "gene_counts": outdir / "gene_counts.tsv",
        "gene_tpm": outdir / "gene_tpm.tsv",
        "transcript_counts": outdir / "transcript_counts.tsv",
        "transcript_tpm": outdir / "transcript_tpm.tsv",
    }
    _write_matrix(written["gene_counts"], "gene_id", samples, gene_counts, integer=True)
    _write_matrix(written["gene_tpm"], "gene_id", samples, gene_tpm, integer=False)
    _write_matrix(
        written["transcript_counts"],
        "transcript_id",
        samples,
        transcript_counts,
        integer=True,
    )
    _write_matrix(
        written["transcript_tpm"],
        "transcript_id",
        samples,
        transcript_tpm,
        integer=False,
    )
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "quant_dirs",
        nargs="+",
        type=Path,
        help="per-sample salmon output directories; the directory name is the sample",
    )
    parser.add_argument("--tx2gene", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    arguments = parser.parse_args(argv)
    try:
        merge(arguments.quant_dirs, arguments.tx2gene, arguments.outdir)
    except MergeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
