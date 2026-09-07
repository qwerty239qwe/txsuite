"""Stream comparison tables using the standard library (also used in stub runs)."""

import argparse
import csv
from pathlib import Path


def collect_de(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    with Path(args.expected).open(encoding="utf-8-sig", newline="") as handle:
        expected = list(csv.DictReader(handle, delimiter="\t"))
    directories = {Path(p).name: Path(p) for p in args.result_dir}
    prefix = [
        "comparison",
        "group_column",
        "group_value",
        "design",
        "reference",
        "test",
        "method",
    ]
    records, sources, fields = [], [], list(prefix)
    for row in expected:
        name, method = row["comparison"], row["method"] or "deseq2"
        row["method"] = method
        directory = directories.get(name)
        result = directory / f"{method}-results.tsv" if directory else None
        significant = directory / "significant-genes.tsv" if directory else None
        success = bool(result and result.is_file() and significant.is_file())
        records.append(
            {
                **row,
                "status": "success" if success else "failed",
                "result": f"de/{name}/{method}-results.tsv" if success else "",
                "significant": f"de/{name}/significant-genes.tsv" if success else "",
                "error": "" if success else "See the Nextflow log for the failed task",
            }
        )
        if success:
            with result.open(encoding="utf-8", newline="") as handle:
                header = csv.DictReader(handle, delimiter="\t").fieldnames or []
            if not header or set(header) & set(prefix):
                raise ValueError(f"Invalid DE table header: {result}")
            fields.extend(field for field in header if field not in fields)
            sources.append((row, result))
    with (output / "comparison-index.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "comparison",
                "status",
                "method",
                "group_column",
                "group_value",
                "design",
                "reference",
                "test",
                "result",
                "significant",
                "error",
            ],
        )
        writer.writeheader()
        writer.writerows(records)
    with (output / "combined-results.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for metadata, result in sources:
            with result.open(encoding="utf-8", newline="") as source:
                for row in csv.DictReader(source, delimiter="\t"):
                    writer.writerow({**row, **{key: metadata[key] for key in prefix}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("expected")
    parser.add_argument("output")
    parser.add_argument("result_dir", nargs="*")
    collect_de(parser.parse_args())
