from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import (
    analysis_command,
    pseudobulk_command,
    validate_samplesheet,
    workflow_command,
)


class Phase2Test(unittest.TestCase):
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
            analysis = analysis_command(
                image="txsuite/single-cell-python:test",
                input_path=matrix,
                outdir=root / "analysis",
                min_genes=1,
                min_cells=1,
                max_mito_pct=20,
                resolution=0.5,
            )
            self.assertEqual(analysis[:3], ["docker", "run", "--rm"])
            self.assertIn("/opt/txsuite/single_cell.py", analysis)
            self.assertIn("analyze", analysis)

            h5ad = root / "data.h5ad"
            h5ad.touch()
            aggregate = pseudobulk_command(
                image="txsuite/single-cell-python:test",
                h5ad=h5ad,
                outdir=root / "de",
                sample_column="sample",
                design="condition",
            )
            self.assertIn("/opt/txsuite/single_cell.py", aggregate)
            self.assertIn("pseudobulk", aggregate)

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
