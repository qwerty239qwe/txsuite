from __future__ import annotations

import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import build_cellranger_image, cellranger_workflow_command


class CellRangerCommandTests(unittest.TestCase):
    def test_build_reference_from_fasta_and_gtf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastqs = root / "fastqs"
            fastqs.mkdir()
            fasta = root / "GRCh38.primary_assembly.genome.fa"
            fasta.touch()
            gtf = root / "gencode.v50.primary_assembly.annotation.gtf"
            gtf.touch()

            command = cellranger_workflow_command(
                DEFAULT_CONFIG,
                genome_name="GRCh38.p14",
                fastqs=fastqs,
                sample="sample1",
                outdir=root / "results",
                fasta=fasta,
                gtf=gtf,
                threads=8,
                memory_gb=32,
                resume=True,
            )

            self.assertEqual(command[:2], ["nextflow", "run"])
            self.assertIn("single_cell_cellranger.nf", command[2])
            self.assertIn("-resume", command)
            self.assertIn("--genome_name", command)
            self.assertIn("GRCh38.p14", command)
            self.assertIn("--fasta", command)
            self.assertIn(str(fasta.resolve()), command)
            self.assertIn("--gtf", command)
            self.assertIn(str(gtf.resolve()), command)
            self.assertIn("--cellranger_threads", command)
            self.assertIn("8", command)
            self.assertNotIn("--cellranger_reference", command)

    def test_reuse_prebuilt_reference(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastqs = root / "fastqs"
            reference = root / "GRCh38.p14"
            fastqs.mkdir()
            reference.mkdir()

            command = cellranger_workflow_command(
                DEFAULT_CONFIG,
                genome_name="GRCh38.p14",
                fastqs=fastqs,
                sample="sample1",
                outdir=root / "results",
                reference=reference,
                cellranger_image="txsuite/cellranger:local",
            )

            self.assertIn("--cellranger_reference", command)
            self.assertIn(str(reference.resolve()), command)
            self.assertNotIn("--fasta", command)
            self.assertNotIn("--gtf", command)
            self.assertIn("--cellranger_image", command)
            self.assertIn("txsuite/cellranger:local", command)

    def test_reference_and_fasta_together_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastqs = root / "fastqs"
            reference = root / "ref"
            fasta = root / "genome.fa"
            gtf = root / "genes.gtf"
            for path in (fastqs, reference):
                path.mkdir()
            fasta.touch()
            gtf.touch()

            with self.assertRaisesRegex(TxSuiteError, "either --reference"):
                cellranger_workflow_command(
                    DEFAULT_CONFIG,
                    genome_name="GRCh38.p14",
                    fastqs=fastqs,
                    sample="sample1",
                    outdir=root / "results",
                    fasta=fasta,
                    gtf=gtf,
                    reference=reference,
                )

    def test_fasta_without_gtf_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastqs = root / "fastqs"
            fasta = root / "genome.fa"
            fastqs.mkdir()
            fasta.touch()

            with self.assertRaisesRegex(TxSuiteError, "both --fasta and --gtf"):
                cellranger_workflow_command(
                    DEFAULT_CONFIG,
                    genome_name="GRCh38.p14",
                    fastqs=fastqs,
                    sample="sample1",
                    outdir=root / "results",
                    fasta=fasta,
                )

    def test_invalid_genome_name_and_sample_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fastqs = root / "fastqs"
            reference = root / "ref"
            for path in (fastqs, reference):
                path.mkdir()

            with self.assertRaisesRegex(TxSuiteError, "genome name"):
                cellranger_workflow_command(
                    DEFAULT_CONFIG,
                    genome_name="bad name!",
                    fastqs=fastqs,
                    sample="sample1",
                    outdir=root / "results",
                    reference=reference,
                )
            with self.assertRaisesRegex(TxSuiteError, "sample"):
                cellranger_workflow_command(
                    DEFAULT_CONFIG,
                    genome_name="GRCh38.p14",
                    fastqs=fastqs,
                    sample="bad sample!",
                    outdir=root / "results",
                    reference=reference,
                )

    def test_missing_inputs_are_rejected_unless_check_inputs_is_false(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing_fastqs = root / "missing-fastqs"
            missing_reference = root / "missing-ref"

            with self.assertRaisesRegex(TxSuiteError, "FASTQ directory"):
                cellranger_workflow_command(
                    DEFAULT_CONFIG,
                    genome_name="GRCh38.p14",
                    fastqs=missing_fastqs,
                    sample="sample1",
                    outdir=root / "results",
                    reference=missing_reference,
                )

            command = cellranger_workflow_command(
                DEFAULT_CONFIG,
                genome_name="GRCh38.p14",
                fastqs=missing_fastqs,
                sample="sample1",
                outdir=root / "results",
                reference=missing_reference,
                check_inputs=False,
            )
            self.assertIn("--cellranger_reference", command)

    def test_build_cellranger_image_requires_an_existing_tarball(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(TxSuiteError, "tarball does not exist"):
                build_cellranger_image(
                    "txsuite/cellranger:test",
                    source_tarball=root / "cellranger-9.0.0.tar.gz",
                    run_dir=root / "run",
                )


class CellRangerPackagedNextflowTests(unittest.TestCase):
    def test_module_and_entrypoint_are_packaged_and_stubbable(self) -> None:
        package = resources.files("txsuite.resources.nextflow")
        module = package.joinpath("modules/cellranger.nf").read_text(encoding="utf-8")
        entrypoint = package.joinpath("single_cell_cellranger.nf").read_text(
            encoding="utf-8"
        )
        config = package.joinpath("nextflow.config").read_text(encoding="utf-8")

        self.assertIn("process CELLRANGER_MKREF", module)
        self.assertIn("process CELLRANGER_COUNT", module)
        self.assertIn("container params.cellranger_image", module)
        self.assertIn("cellranger mkref", module)
        self.assertIn("cellranger count", module)
        self.assertEqual(module.count("stub:"), 2)

        self.assertIn("CELLRANGER_MKREF(fasta_ch, gtf_ch)", entrypoint)
        self.assertIn("CELLRANGER_COUNT(reference_ch, fastqs_ch)", entrypoint)

        for name in ("genome_name", "cellranger_reference", "fastqs", "sample", "cellranger_image"):
            self.assertIn(name, config)


if __name__ == "__main__":
    unittest.main()
