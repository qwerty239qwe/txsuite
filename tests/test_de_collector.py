import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from txsuite.resources.single_cell_python.collect_de import collect_de


class CollectorTest(unittest.TestCase):
    def test_streams_mixed_tables_and_reports_missing_comparisons(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expected = root / "expected.tsv"
            expected.write_text(
                "comparison\tgroup_column\tgroup_value\tdesign\treference\ttest\tmethod\n"
                "a\t\t\tcondition\tcontrol\ttreated\tdeseq2\n"
                "b\t\t\tcondition\tcontrol\ttreated\tedger\n"
                "c\t\t\tcondition\tcontrol\ttreated\tlimma\n",
                encoding="utf-8",
            )
            for name, method, column in (
                ("a", "deseq2", "baseMean"),
                ("b", "edger", "logCPM"),
            ):
                directory = root / name
                directory.mkdir()
                (directory / f"{method}-results.tsv").write_text(
                    f"gene_id\t{column}\tpadj\n001\t1.20\tNA\n", encoding="utf-8"
                )
                (directory / "significant-genes.tsv").write_text(
                    "gene_id\n", encoding="utf-8"
                )
            for directories in ([root / "a", root / "b"], []):
                with self.subTest(directories=directories):
                    collect_de(
                        SimpleNamespace(
                            expected=expected,
                            output=root / "out",
                            result_dir=directories,
                        )
                    )
                    with (root / "out/comparison-index.tsv").open(newline="") as handle:
                        index = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual(
                        [r["status"] for r in index],
                        ["success", "success", "failed"]
                        if directories
                        else ["failed"] * 3,
                    )
                    with (root / "out/combined-results.tsv").open(newline="") as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual(len(rows), len(directories))
                    if rows:
                        self.assertEqual(rows[0]["gene_id"], "001")
                        self.assertEqual(rows[0]["padj"], "NA")
                        self.assertEqual(rows[0]["baseMean"], "1.20")
                        self.assertEqual(rows[0]["logCPM"], "")
