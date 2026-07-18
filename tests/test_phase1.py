from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from txsuite.bulk import deseq2_command, validate_samplesheet, workflow_command
from txsuite.config import DEFAULT_CONFIG
from txsuite.runtime import TxSuiteError, run_command


class Phase1Test(unittest.TestCase):
    def test_bulk_commands_and_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            samplesheet = root / "samples.csv"
            samplesheet.write_text(
                "sample,fastq_1,fastq_2,strandedness\n"
                "control,control_R1.fastq.gz,control_R2.fastq.gz,auto\n",
                encoding="utf-8",
            )
            self.assertEqual(validate_samplesheet(samplesheet), 1)
            command = workflow_command(
                DEFAULT_CONFIG,
                samplesheet=samplesheet,
                outdir=root / "results",
                resume=True,
            )
            self.assertEqual(
                command[:7],
                [
                    "nextflow",
                    "run",
                    "nf-core/rnaseq",
                    "-r",
                    "3.26.0",
                    "-profile",
                    "docker",
                ],
            )
            self.assertEqual(command[-1], "-resume")

            counts = root / "counts.tsv"
            metadata = root / "metadata.tsv"
            counts.write_text("gene_id\tcontrol\nA\t1\n", encoding="utf-8")
            metadata.write_text(
                "sample\tcondition\ncontrol\tuntreated\n", encoding="utf-8"
            )
            de = deseq2_command(
                image="txsuite/bulk-r:test",
                counts=counts,
                metadata=metadata,
                design="condition",
                reference="untreated",
                test="treated",
                outdir=root / "de",
            )
            self.assertEqual(de[:3], ["docker", "run", "--rm"])

            run_dir = root / "run"
            run_command(
                [sys.executable, "-c", "print('ok')"],
                run_dir=run_dir,
                task="test",
                backend="python",
                inputs={},
                outputs={},
                artifacts=[],
            )
            self.assertEqual(
                json.loads((run_dir / "run.json").read_text())["status"], "success"
            )
            self.assertTrue((run_dir / "txsuite-results.json").exists())

    def test_samplesheet_rejects_missing_required_column(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.csv"
            path.write_text("sample,fastq_1\na,a.fastq.gz\n", encoding="utf-8")
            with self.assertRaisesRegex(TxSuiteError, "missing columns"):
                validate_samplesheet(path)


if __name__ == "__main__":
    unittest.main()
