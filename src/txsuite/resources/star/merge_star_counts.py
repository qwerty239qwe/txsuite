"""Merge STAR per-sample gene counts, inferring library strandedness.

``STAR --quantMode GeneCounts`` writes ``ReadsPerGene.out.tab`` with one row per
gene and three count columns: unstranded, forward, and reverse. Choosing the
wrong column does not fail — it silently produces a matrix that is mostly noise,
which is why the choice is inferred from the data, recorded alongside the
matrix, and refused rather than guessed when the evidence is weak.

The first four rows are alignment summary counters rather than genes. They are
excluded from both the matrix and the inference totals; leaving them in would
let ``N_noFeature`` dominate the column sums and invert the answer.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

COUNTS_FILE = "ReadsPerGene.out.tab"
SUMMARY_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
# Column order in ReadsPerGene.out.tab after the gene ID.
STRAND_COLUMNS = ("unstranded", "forward", "reverse")
# A stranded library sends at least this fraction of gene reads to one
# orientation; an unstranded library lands within this band of an even split.
DEFAULT_STRANDED_FRACTION = 0.8
DEFAULT_UNSTRANDED_BAND = 0.15


class MergeError(RuntimeError):
    """A STAR count file is missing, malformed, or its strandedness is unclear."""


def read_star_counts(path: Path) -> tuple[dict[str, tuple[int, int, int]], dict[str, int]]:
    """Return per-gene counts and the excluded alignment summary counters."""

    genes: dict[str, tuple[int, int, int]] = {}
    summary: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
            if not row or not row[0].strip():
                continue
            if len(row) < 4:
                raise MergeError(
                    f"{path}:{line_number}: expected a gene ID and three count columns"
                )
            name = row[0].strip()
            try:
                values = (int(row[1]), int(row[2]), int(row[3]))
            except (TypeError, ValueError) as exc:
                raise MergeError(f"{path}:{line_number}: non-integer count") from exc
            if name in SUMMARY_ROWS:
                summary[name] = values[0]
                continue
            if name in genes:
                raise MergeError(f"{path}:{line_number}: repeated gene {name!r}")
            genes[name] = values
    if not genes:
        raise MergeError(f"No gene rows in {path}")
    return genes, summary


def infer_strandedness(
    genes: dict[str, tuple[int, int, int]],
    *,
    stranded_fraction: float = DEFAULT_STRANDED_FRACTION,
    unstranded_band: float = DEFAULT_UNSTRANDED_BAND,
) -> tuple[str, tuple[int, int, int], float]:
    """Classify the library from the balance between the two stranded columns.

    The unstranded column is not a competitor: STAR fills it with reads counted
    regardless of orientation, so it is approximately forward + reverse and
    would always win a naive "largest column" comparison. The signal is the
    *split* between forward and reverse. A stranded library sends nearly all
    reads to one of them; an unstranded library splits them near evenly.
    """

    totals = tuple(sum(values[index] for values in genes.values()) for index in range(3))
    forward, reverse = totals[1], totals[2]
    stranded_total = forward + reverse
    if stranded_total == 0:
        raise MergeError("no reads assigned to genes in either stranded column")

    forward_fraction = forward / stranded_total
    if forward_fraction >= stranded_fraction:
        return "forward", totals, forward_fraction
    if forward_fraction <= 1.0 - stranded_fraction:
        return "reverse", totals, forward_fraction
    if abs(forward_fraction - 0.5) <= unstranded_band:
        return "unstranded", totals, forward_fraction
    raise MergeError(
        "cannot infer strandedness: forward/reverse split is "
        f"{forward_fraction:.2f} (forward={forward}, reverse={reverse}), which is "
        "neither clearly stranded nor evenly split; pass an explicit strandedness"
    )


def merge(
    quant_dirs: list[Path],
    outdir: Path,
    *,
    strandedness: str = "auto",
    stranded_fraction: float = DEFAULT_STRANDED_FRACTION,
) -> dict[str, Path]:
    """Write the gene-count matrix and the strandedness report."""

    if strandedness not in ("auto",) + STRAND_COLUMNS:
        raise MergeError(f"unknown strandedness: {strandedness!r}")
    if not quant_dirs:
        raise MergeError("At least one STAR output directory is required")

    samples: list[str] = []
    per_sample: dict[str, dict[str, tuple[int, int, int]]] = {}
    report: list[dict[str, object]] = []
    inferred: set[str] = set()

    for quant_dir in quant_dirs:
        # Accepts either a per-sample directory or STAR's prefixed file, since
        # the DAG names outputs "<sample>.ReadsPerGene.out.tab" so that they stay
        # unique at the top level of a task directory.
        if quant_dir.is_dir():
            sample, path = quant_dir.name, quant_dir / COUNTS_FILE
        else:
            sample = quant_dir.name
            if sample.endswith(f".{COUNTS_FILE}"):
                sample = sample[: -len(COUNTS_FILE) - 1]
            path = quant_dir
        if sample in per_sample:
            raise MergeError(f"Duplicate sample name: {sample}")
        if not path.is_file():
            raise MergeError(f"No {COUNTS_FILE} under {quant_dir}")
        genes, summary = read_star_counts(path)
        try:
            choice, totals, fraction = infer_strandedness(
                genes, stranded_fraction=stranded_fraction
            )
        except MergeError as exc:
            raise MergeError(f"{sample}: {exc}") from exc
        if strandedness != "auto":
            choice = strandedness
        inferred.add(choice)
        samples.append(sample)
        per_sample[sample] = genes
        report.append(
            {
                "sample": sample,
                "strandedness": choice,
                "source": "inferred" if strandedness == "auto" else "explicit",
                "unstranded": totals[0],
                "forward": totals[1],
                "reverse": totals[2],
                "forward_fraction": f"{fraction:.4f}",
                "unmapped": summary.get("N_unmapped", 0),
                "no_feature": summary.get("N_noFeature", 0),
            }
        )

    if len(inferred) > 1:
        raise MergeError(
            "samples disagree on strandedness ("
            + ", ".join(f"{row['sample']}={row['strandedness']}" for row in report)
            + "); mixed library preparations cannot share one count matrix"
        )
    column = STRAND_COLUMNS.index(next(iter(inferred)))

    outdir.mkdir(parents=True, exist_ok=True)
    written = {
        "counts": outdir / "gene_counts.tsv",
        "strandedness": outdir / "strandedness.tsv",
    }
    gene_ids = sorted({gene for genes in per_sample.values() for gene in genes})
    with written["counts"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", *samples])
        for gene in gene_ids:
            writer.writerow(
                [gene]
                + [str(per_sample[s].get(gene, (0, 0, 0))[column]) for s in samples]
            )
    with written["strandedness"].open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "sample",
            "strandedness",
            "source",
            "unstranded",
            "forward",
            "reverse",
            "forward_fraction",
            "unmapped",
            "no_feature",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(report)
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "quant_dirs",
        nargs="+",
        type=Path,
        help="per-sample STAR output directories; the directory name is the sample",
    )
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--strandedness", default="auto", choices=("auto",) + STRAND_COLUMNS
    )
    parser.add_argument(
        "--stranded-fraction", type=float, default=DEFAULT_STRANDED_FRACTION
    )
    arguments = parser.parse_args(argv)
    try:
        merge(
            arguments.quant_dirs,
            arguments.outdir,
            strandedness=arguments.strandedness,
            stranded_fraction=arguments.stranded_fraction,
        )
    except MergeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
