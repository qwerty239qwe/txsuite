from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

from txsuite.bulk import workflow_command
from txsuite.config import DEFAULT_CONFIG
from txsuite.hardening import cache_reference, image_is_locked
from txsuite.runtime import TxSuiteError


class Phase4Test(unittest.TestCase):
    def test_verified_cache_digest_gate_and_nextflow_config(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "reference.fa.gz"
            source.write_bytes(b"reference-data")
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            cached = cache_reference(
                str(source), digest, "GRCh38.fa.gz", root / "cache"
            )
            self.assertEqual(cached.read_bytes(), b"reference-data")
            self.assertTrue((root / "cache" / "GRCh38.fa.gz.json").is_file())
            with self.assertRaisesRegex(TxSuiteError, "checksum mismatch"):
                cache_reference(str(source), "0" * 64, "bad.fa.gz", root / "cache")
            self.assertFalse((root / "cache" / "bad.fa.gz").exists())

            self.assertTrue(image_is_locked("ghcr.io/org/image@sha256:" + "a" * 64))
            self.assertTrue(image_is_locked("sha256:" + "b" * 64))
            self.assertFalse(image_is_locked("txsuite/image:0.1.0"))

            samplesheet = root / "samples.csv"
            samplesheet.write_text(
                "sample,fastq_1,fastq_2,strandedness\n"
                "sample,a_R1.fastq.gz,a_R2.fastq.gz,auto\n",
                encoding="utf-8",
            )
            nextflow_config = root / "slurm.config"
            nextflow_config.write_text("process.executor = 'slurm'\n", encoding="utf-8")
            command = workflow_command(
                DEFAULT_CONFIG,
                samplesheet=samplesheet,
                outdir=root / "results",
                nextflow_config=nextflow_config,
                resume=True,
            )
            self.assertEqual(
                command[-3:], ["-c", str(nextflow_config.resolve()), "-resume"]
            )


if __name__ == "__main__":
    unittest.main()
