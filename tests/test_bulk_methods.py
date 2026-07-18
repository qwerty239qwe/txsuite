from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
import unittest

from txsuite.bulk import deseq2_command, enrichment_command
from txsuite.cli import run
from txsuite.runtime import TxSuiteError


class BulkMethodTests(unittest.TestCase):
    def test_de_and_enrichment_commands_expose_the_advanced_methods(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts = root / "counts.tsv"
            metadata = root / "metadata.tsv"
            de_results = root / "de.tsv"
            genesets = root / "sets.gmt"
            counts.write_text("gene_id\ta\tb\nG1\t1\t2\n", encoding="utf-8")
            metadata.write_text(
                "sample\tcondition\tbatch\na\tcontrol\t1\nb\ttreated\t2\n",
                encoding="utf-8",
            )
            de_results.write_text(
                "gene_id\tlog2FoldChange\tstat\tpadj\nG1\t2\t4\t0.01\n",
                encoding="utf-8",
            )
            genesets.write_text("SET\tdescription\tG1\n", encoding="utf-8")

            de = deseq2_command(
                image="txsuite/bulk-r:test",
                counts=counts,
                metadata=metadata,
                design="condition",
                reference="control",
                test="treated",
                outdir=root / "de",
                covariates=("batch",),
                padj=0.01,
                lfc=1.5,
                top_genes=20,
            )
            self.assertIn("/opt/txsuite/deseq2.R", de)
            self.assertEqual(de[-4:], ["0.01", "1.5", "20", "batch"])

            enrich = enrichment_command(
                image="txsuite/bulk-r:test",
                de_results=de_results,
                genesets=genesets,
                mode="gsea",
                outdir=root / "enrichment",
                min_size=1,
            )
            self.assertIn("/opt/txsuite/enrichment.R", enrich)
            self.assertIn("gsea", enrich)

            with self.assertRaisesRegex(TxSuiteError, "different from design"):
                deseq2_command(
                    image="txsuite/bulk-r:test",
                    counts=counts,
                    metadata=metadata,
                    design="condition",
                    reference="control",
                    test="treated",
                    outdir=root / "de",
                    covariates=("condition",),
                )

            output = StringIO()
            with redirect_stdout(output):
                status = run(
                    [
                        "bulk",
                        "enrich",
                        "--de",
                        str(de_results),
                        "--genesets",
                        str(genesets),
                        "--mode",
                        "ora",
                        "--min-size",
                        "1",
                        "--outdir",
                        str(root / "ora"),
                        "--dry-run",
                    ]
                )
            self.assertEqual(status, 0)
            self.assertIn("enrichment.R", output.getvalue())

        resource_root = (
            Path(__file__).parents[1] / "src" / "txsuite" / "resources" / "bulk_r"
        )
        deseq2 = (resource_root / "deseq2.R").read_text(encoding="utf-8")
        enrichment = (resource_root / "enrichment.R").read_text(encoding="utf-8")
        for output_name in (
            "sample-qc.tsv",
            "pca.pdf",
            "ma.pdf",
            "volcano.pdf",
            "top-genes-heatmap.pdf",
        ):
            self.assertIn(output_name, deseq2)
        self.assertIn("enricher(", enrichment)
        self.assertIn("GSEA(", enrichment)


if __name__ == "__main__":
    unittest.main()
