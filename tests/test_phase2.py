from __future__ import annotations

import json
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from txsuite.cli import run
from txsuite.config import DEFAULT_CONFIG
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import (
    analysis_command,
    pseudobulk_command,
    run_pseudobulk_manifest,
    validate_samplesheet,
    workflow_command,
)


class Phase2Test(unittest.TestCase):
    def test_manifest_pseudobulk_runs_combines_and_resumes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "data.h5ad"
            h5ad.touch()
            manifest = root / "comparisons.tsv"
            manifest.write_text(
                "comparison\tgroup_column\tgroup_value\tdesign\treference\ttest\n"
                "t_cells\tcell_type\tT_cell\tcondition\tcontrol\ttreated\n",
                encoding="utf-8",
            )
            outdir = root / "batch"

            def fake_run(command, **kwargs) -> None:
                run_dir = kwargs["run_dir"]
                run_dir.mkdir(parents=True, exist_ok=True)
                (run_dir / "run.json").write_text(
                    json.dumps({"status": "success"}), encoding="utf-8"
                )
                outputs = kwargs["outputs"]
                if "counts" in outputs:
                    Path(outputs["counts"]).write_text(
                        "gene_id\ta\tb\nG1\t1\t2\n", encoding="utf-8"
                    )
                    Path(outputs["metadata"]).write_text(
                        "sample\tcondition\na\tcontrol\nb\ttreated\n",
                        encoding="utf-8",
                    )
                else:
                    result_dir = Path(outputs["outdir"])
                    table = (
                        "gene_id\tbaseMean\tlog2FoldChange\tpadj\nG1\t1.5\t1\t0.05\n"
                    )
                    (result_dir / "deseq2-results.tsv").write_text(
                        table, encoding="utf-8"
                    )
                    (result_dir / "significant-genes.tsv").write_text(
                        table, encoding="utf-8"
                    )

            with patch("txsuite.single_cell.run_command", side_effect=fake_run):
                run_pseudobulk_manifest(
                    manifest=manifest,
                    h5ad=h5ad,
                    outdir=outdir,
                    sample_column="sample",
                    single_cell_image="txsuite/single-cell:test",
                    bulk_image="txsuite/bulk:test",
                )
            self.assertIn("t_cells", (outdir / "combined-results.tsv").read_text())
            self.assertIn("success", (outdir / "comparison-index.tsv").read_text())

            with patch("txsuite.single_cell.run_command") as execute:
                run_pseudobulk_manifest(
                    manifest=manifest,
                    h5ad=h5ad,
                    outdir=outdir,
                    sample_column="sample",
                    single_cell_image="txsuite/single-cell:test",
                    bulk_image="txsuite/bulk:test",
                    resume=True,
                )
            execute.assert_not_called()
            self.assertIn("skipped", (outdir / "comparison-index.tsv").read_text())

            output = StringIO()
            with redirect_stdout(output):
                status = run(
                    [
                        "single-cell",
                        "pseudobulk-batch",
                        "--input",
                        str(h5ad),
                        "--sample-column",
                        "sample",
                        "--manifest",
                        str(manifest),
                        "--outdir",
                        str(outdir),
                        "--image",
                        "txsuite/single-cell:test",
                        "--bulk-image",
                        "txsuite/bulk:test",
                        "--direct",
                        "--dry-run",
                    ]
                )
            self.assertEqual(status, 0)
            self.assertEqual(output.getvalue().count("docker run"), 2)

            output = StringIO()
            with redirect_stdout(output):
                status = run(
                    [
                        "single-cell",
                        "pseudobulk-batch",
                        "--input",
                        str(h5ad),
                        "--sample-column",
                        "sample",
                        "--manifest",
                        str(manifest),
                        "--outdir",
                        str(outdir),
                        "--dry-run",
                    ]
                )
            self.assertEqual(status, 0)
            self.assertIn("nextflow run", output.getvalue())
            self.assertIn("--manifest", output.getvalue())

    def test_single_cell_workflow_and_downstream_commands(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            samplesheet = root / "samples.csv"
            samplesheet.write_text(
                "sample,fastq_1,fastq_2\n"
                "control,control_R1.fastq.gz,control_R2.fastq.gz\n",
                encoding="utf-8",
            )
            self.assertEqual(validate_samplesheet(samplesheet), 1)
            workflow = workflow_command(
                DEFAULT_CONFIG,
                samplesheet=samplesheet,
                outdir=root / "workflow",
                aligner="star",
                protocol="10XV3",
                resume=True,
            )
            self.assertEqual(
                workflow[:7],
                [
                    "nextflow",
                    "run",
                    "nf-core/scrnaseq",
                    "-r",
                    "4.2.0",
                    "-profile",
                    "docker",
                ],
            )
            self.assertIn("star", workflow)
            self.assertEqual(workflow[-1], "-resume")

            matrix = root / "matrix"
            matrix.mkdir()
            cell_metadata = root / "metadata.tsv"
            cell_metadata.write_text(
                "barcode\tsample\tcondition\tbatch\n"
                "cell-1\tsample-1\tcontrol\tbatch-1\n",
                encoding="utf-8",
            )
            analysis = analysis_command(
                image="txsuite/single-cell-python:test",
                input_path=matrix,
                outdir=root / "analysis",
                min_genes=1,
                min_cells=1,
                max_mito_pct=20,
                resolution=0.5,
                metadata=cell_metadata,
                batch_column="batch",
                integration="harmony",
                doublets="score",
                doublet_batch_column="sample",
                expected_doublet_rate=0.1,
                doublet_threshold=0.2,
                top_markers=25,
                counts_layer="raw_counts",
                target_sum=5_000,
                n_hvg=1_000,
                hvg_flavor="cell_ranger",
                n_pcs=30,
                n_neighbors=10,
                umap_min_dist=0.2,
                marker_method="t-test",
            )
            self.assertEqual(analysis[:3], ["docker", "run", "--rm"])
            self.assertIn("/opt/txsuite/single_cell.py", analysis)
            self.assertIn("analyze", analysis)
            self.assertIn("/input/metadata.tsv", analysis)
            self.assertIn("--batch-column", analysis)
            self.assertIn("--integration", analysis)
            self.assertEqual(analysis[analysis.index("--integration") + 1], "harmony")
            self.assertIn("--doublet-batch-column", analysis)
            self.assertIn("--top-markers", analysis)
            self.assertEqual(
                analysis[analysis.index("--counts-layer") + 1], "raw_counts"
            )
            self.assertEqual(
                analysis[analysis.index("--hvg-flavor") + 1], "cell_ranger"
            )
            self.assertEqual(analysis[analysis.index("--marker-method") + 1], "t-test")

            with self.assertRaisesRegex(TxSuiteError, "requires a batch column"):
                analysis_command(
                    image="txsuite/single-cell-python:test",
                    input_path=matrix,
                    outdir=root / "analysis",
                    min_genes=1,
                    min_cells=1,
                    max_mito_pct=20,
                    resolution=0.5,
                    integration="harmony",
                )
            with self.assertRaisesRegex(TxSuiteError, "must be none or harmony"):
                analysis_command(
                    image="txsuite/single-cell-python:test",
                    input_path=matrix,
                    outdir=root / "analysis",
                    min_genes=1,
                    min_cells=1,
                    max_mito_pct=20,
                    resolution=0.5,
                    integration="unknown",
                )
            with self.assertRaisesRegex(TxSuiteError, "requires PCA"):
                analysis_command(
                    image="txsuite/single-cell-python:test",
                    input_path=matrix,
                    outdir=root / "analysis",
                    min_genes=1,
                    min_cells=1,
                    max_mito_pct=20,
                    resolution=0.5,
                    batch_column="batch",
                    integration="harmony",
                    stop_after="qc",
                )

            h5ad = root / "data.h5ad"
            h5ad.touch()
            aggregate = pseudobulk_command(
                image="txsuite/single-cell-python:test",
                h5ad=h5ad,
                outdir=root / "de",
                sample_column="sample",
                design="condition",
                group_column="cell_type",
                group_value="T_cell",
                covariates=("batch",),
                reference="control",
                test="treated",
            )
            self.assertIn("/opt/txsuite/single_cell.py", aggregate)
            self.assertIn("pseudobulk", aggregate)
            self.assertIn("--group-column", aggregate)
            self.assertIn("--covariate", aggregate)

            with self.assertRaisesRegex(TxSuiteError, "used together"):
                pseudobulk_command(
                    image="txsuite/single-cell-python:test",
                    h5ad=h5ad,
                    outdir=root / "de",
                    sample_column="sample",
                    design="condition",
                    group_column="cell_type",
                )

    def test_samplesheet_requires_official_column_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.csv"
            path.write_text(
                "fastq_1,sample,fastq_2\na.fastq.gz,a,a_R2.fastq.gz\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(TxSuiteError, "first columns"):
                validate_samplesheet(path)


if __name__ == "__main__":
    unittest.main()
