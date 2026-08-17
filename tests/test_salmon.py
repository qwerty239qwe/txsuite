from __future__ import annotations

import importlib.util
import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.bulk import salmon_workflow_command, workflow_command
from txsuite.config import load_config
from txsuite.project.adapters.bulk import bulk_salmon_command
from txsuite.project.config import load_project_config
from txsuite.project.planner import plan_workflow
from txsuite.project.presets import scaffold_project_preset
from txsuite.project.adapters.single_cell import alevin_command
from txsuite.project.registry import get_stage_spec, validate_stage_parameters
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import alevin_workflow_command
from txsuite.single_cell import workflow_command as scrnaseq_workflow_command


def _load_resource_module(package: str, filename: str, name: str):
    path = Path(str(resources.files(package).joinpath(filename)))
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SalmonWorkflowCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.samplesheet = self.root / "samplesheet.csv"
        self.samplesheet.write_text(
            "sample,fastq_1,fastq_2\nsample1,r1.fastq.gz,r2.fastq.gz\n", encoding="utf-8"
        )
        self.fasta = self.root / "genome.fa"
        self.fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.gtf = self.root / "genes.gtf"
        self.gtf.write_text("#gtf\n", encoding="utf-8")
        self.index = self.root / "salmon_index"
        self.index.mkdir()
        self.tx2gene = self.root / "tx2gene.tsv"
        self.tx2gene.write_text("tx1\tgene1\n", encoding="utf-8")
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def test_reference_build_command_carries_pinned_image_and_workflow(self) -> None:
        command = salmon_workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            fasta=self.fasta,
            gtf=self.gtf,
        )

        self.assertEqual(command[:3], ["nextflow", "run", command[2]])
        self.assertTrue(command[2].endswith("bulk_salmon.nf"))
        self.assertEqual(command[command.index("-profile") + 1], "docker")
        self.assertEqual(
            command[command.index("--samplesheet") + 1], str(self.samplesheet.resolve())
        )
        self.assertEqual(command[command.index("--fasta") + 1], str(self.fasta.resolve()))
        self.assertEqual(command[command.index("--gtf") + 1], str(self.gtf.resolve()))
        self.assertEqual(
            command[command.index("--salmon_image") + 1], self.config["images"]["salmon"]
        )
        self.assertEqual(command[command.index("--salmon_libtype") + 1], "A")
        self.assertNotIn("--salmon_index", command)
        self.assertNotIn("-resume", command)

    def test_prebuilt_index_requires_tx2gene_and_omits_reference_inputs(self) -> None:
        command = salmon_workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            salmon_index=self.index,
            tx2gene=self.tx2gene,
            libtype="ISR",
            resume=True,
        )

        self.assertEqual(
            command[command.index("--salmon_index") + 1], str(self.index.resolve())
        )
        self.assertEqual(command[command.index("--tx2gene") + 1], str(self.tx2gene.resolve()))
        self.assertEqual(command[command.index("--salmon_libtype") + 1], "ISR")
        self.assertNotIn("--fasta", command)
        self.assertIn("-resume", command)

        with self.assertRaises(TxSuiteError):
            salmon_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                salmon_index=self.index,
            )

    def test_reference_selection_is_exclusive_and_complete(self) -> None:
        for kwargs in (
            {},
            {"fasta": None, "gtf": None},
            {"fasta": None, "gtf": None, "salmon_index": None},
            {"fasta": None},
            {"salmon_index": None},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                salmon_workflow_command(
                    self.config,
                    samplesheet=self.samplesheet,
                    outdir=self.outdir,
                    **kwargs,
                )

        with self.assertRaises(TxSuiteError):
            salmon_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                fasta=self.fasta,
                gtf=self.gtf,
                salmon_index=self.index,
                tx2gene=self.tx2gene,
            )

        with self.assertRaises(TxSuiteError):
            salmon_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                fasta=self.fasta,
            )

    def test_missing_inputs_and_bad_parameters_are_rejected(self) -> None:
        with self.assertRaises(TxSuiteError):
            salmon_workflow_command(
                self.config,
                samplesheet=self.root / "absent.csv",
                outdir=self.outdir,
                fasta=self.fasta,
                gtf=self.gtf,
            )
        with self.assertRaises(TxSuiteError):
            salmon_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                fasta=self.fasta,
                gtf=self.gtf,
                libtype="not-a-libtype",
            )
        with self.assertRaises(TxSuiteError):
            salmon_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                fasta=self.fasta,
                gtf=self.gtf,
                kmer_len=0,
            )


class AlevinWorkflowCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.samplesheet = self.root / "samplesheet.csv"
        self.samplesheet.write_text(
            "sample,fastq_1,fastq_2\nsample1,r1.fastq.gz,r2.fastq.gz\n", encoding="utf-8"
        )
        self.fasta = self.root / "genome.fa"
        self.fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.gtf = self.root / "genes.gtf"
        self.gtf.write_text("#gtf\n", encoding="utf-8")
        self.index = self.root / "simpleaf_index"
        self.index.mkdir()
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def test_command_names_both_images_and_the_native_dag(self) -> None:
        command = alevin_workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            fasta=self.fasta,
            gtf=self.gtf,
        )

        self.assertTrue(command[2].endswith("single_cell_alevin.nf"))
        self.assertEqual(
            command[command.index("--salmon_image") + 1], self.config["images"]["salmon"]
        )
        self.assertEqual(
            command[command.index("--single_cell_image") + 1],
            self.config["images"]["single_cell_python"],
        )
        self.assertEqual(command[command.index("--alevin_chemistry") + 1], "10xv3")
        self.assertEqual(command[command.index("--alevin_resolution") + 1], "cr-like")
        self.assertNotIn("--alevin_whitelist", command)

    def test_index_and_reference_are_mutually_exclusive(self) -> None:
        command = alevin_workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            simpleaf_index=self.index,
            chemistry="10xv2",
        )
        self.assertEqual(
            command[command.index("--simpleaf_index") + 1], str(self.index.resolve())
        )
        self.assertEqual(command[command.index("--alevin_chemistry") + 1], "10xv2")

        with self.assertRaises(TxSuiteError):
            alevin_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                simpleaf_index=self.index,
                fasta=self.fasta,
                gtf=self.gtf,
            )
        with self.assertRaises(TxSuiteError):
            alevin_workflow_command(
                self.config, samplesheet=self.samplesheet, outdir=self.outdir
            )

    def test_unknown_chemistry_and_resolution_are_rejected(self) -> None:
        for kwargs in (
            {"chemistry": "10xv99"},
            {"resolution": "guess"},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                alevin_workflow_command(
                    self.config,
                    samplesheet=self.samplesheet,
                    outdir=self.outdir,
                    simpleaf_index=self.index,
                    **kwargs,
                )

    def test_whitelist_is_passed_through_when_present(self) -> None:
        whitelist = self.root / "whitelist.txt"
        whitelist.write_text("AAACCTGAGAAACCAT\n", encoding="utf-8")
        command = alevin_workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            simpleaf_index=self.index,
            whitelist=whitelist,
        )
        self.assertEqual(
            command[command.index("--alevin_whitelist") + 1], str(whitelist.resolve())
        )


class NfCorePassThroughTests(unittest.TestCase):
    """The pinned nf-core pipelines already run salmon; these expose its knobs."""

    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.samplesheet = self.root / "samplesheet.csv"
        self.samplesheet.write_text("sample\n", encoding="utf-8")
        self.index = self.root / "index"
        self.index.mkdir()
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def test_bulk_rnaseq_forwards_pseudo_aligner_settings(self) -> None:
        command = workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            pseudo_aligner="salmon",
            skip_alignment=True,
            salmon_index=self.index,
        )
        self.assertEqual(command[command.index("--pseudo_aligner") + 1], "salmon")
        self.assertIn("--skip_alignment", command)
        self.assertEqual(
            command[command.index("--salmon_index") + 1], str(self.index.resolve())
        )

        baseline = workflow_command(
            self.config, samplesheet=self.samplesheet, outdir=self.outdir
        )
        self.assertNotIn("--pseudo_aligner", baseline)
        self.assertNotIn("--skip_alignment", baseline)

    def test_bulk_rnaseq_rejects_inconsistent_pseudo_aligner_settings(self) -> None:
        for kwargs in (
            {"pseudo_aligner": "bowtie"},
            {"skip_alignment": True},
            {"salmon_index": self.index},
            {"pseudo_aligner": "kallisto", "salmon_index": self.index},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                workflow_command(
                    self.config,
                    samplesheet=self.samplesheet,
                    outdir=self.outdir,
                    **kwargs,
                )

    def test_scrnaseq_forwards_simpleaf_index_and_txp2gene(self) -> None:
        txp2gene = self.root / "txp2gene.tsv"
        txp2gene.write_text("tx1\tgene1\n", encoding="utf-8")
        command = scrnaseq_workflow_command(
            self.config,
            samplesheet=self.samplesheet,
            outdir=self.outdir,
            simpleaf_index=self.index,
            txp2gene=txp2gene,
        )
        self.assertEqual(
            command[command.index("--simpleaf_index") + 1], str(self.index.resolve())
        )
        self.assertEqual(command[command.index("--txp2gene") + 1], str(txp2gene.resolve()))

        with self.assertRaises(TxSuiteError):
            scrnaseq_workflow_command(
                self.config,
                samplesheet=self.samplesheet,
                outdir=self.outdir,
                aligner="star",
                simpleaf_index=self.index,
            )

    def test_registry_defaults_keep_the_upstream_behaviour_unchanged(self) -> None:
        bulk_defaults = validate_stage_parameters(get_stage_spec("bulk.rnaseq"), {})
        self.assertIsNone(bulk_defaults["pseudo_aligner"])
        self.assertFalse(bulk_defaults["skip_alignment"])
        single_cell_defaults = validate_stage_parameters(
            get_stage_spec("single-cell.scrnaseq"), {}
        )
        self.assertIsNone(single_cell_defaults["simpleaf_index"])
        with self.assertRaises(TxSuiteError):
            validate_stage_parameters(
                get_stage_spec("bulk.rnaseq"), {"pseudo_aligner": "bowtie"}
            )


class SalmonStageRegistryTests(unittest.TestCase):
    def test_bulk_salmon_is_a_drop_in_producer_for_bulk_de(self) -> None:
        spec = get_stage_spec("bulk.salmon")
        de_spec = get_stage_spec("bulk.de")

        self.assertEqual(spec.modality, "bulk")
        self.assertEqual(spec.inputs, get_stage_spec("bulk.rnaseq").inputs)
        self.assertEqual(spec.outputs["counts"], de_spec.inputs["counts"])
        self.assertEqual(spec.outputs["tx_counts"], "bulk.transcript-counts")
        self.assertEqual(spec.required_images, ("images.salmon",))
        self.assertEqual(spec.required_executables, ("nextflow",))
        self.assertTrue(spec.supports_resume)

        defaults = validate_stage_parameters(spec, {})
        self.assertEqual(defaults["libtype"], "A")
        self.assertEqual(defaults["kmer_len"], 31)
        self.assertIs(defaults["salmon_index"], None)
        with self.assertRaises(TxSuiteError):
            validate_stage_parameters(spec, {"libtype": "nonsense"})

    def test_alevin_is_a_drop_in_producer_for_scanpy(self) -> None:
        spec = get_stage_spec("single-cell.alevin")

        self.assertEqual(spec.modality, "single-cell")
        self.assertEqual(spec.inputs, get_stage_spec("single-cell.scrnaseq").inputs)
        self.assertEqual(spec.outputs["matrix"], get_stage_spec("single-cell.scanpy").inputs["input"])
        self.assertEqual(spec.required_images, ("images.salmon", "images.single_cell_python"))
        self.assertTrue(spec.supports_resume)

        defaults = validate_stage_parameters(spec, {})
        self.assertEqual(defaults["chemistry"], "10xv3")
        with self.assertRaises(TxSuiteError):
            validate_stage_parameters(spec, {"chemistry": "10xv99"})

    def test_adapters_build_the_same_argv_as_the_library_helpers(self) -> None:
        config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            samplesheet = root / "samplesheet.csv"
            samplesheet.write_text(
                "sample,fastq_1,fastq_2\nsample1,r1.fastq.gz,r2.fastq.gz\n",
                encoding="utf-8",
            )
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
            gtf = root / "genes.gtf"
            gtf.write_text("#gtf\n", encoding="utf-8")
            outdir = root / "results"

            context = {
                "global_config": config,
                "resolved_inputs": {"samplesheet": samplesheet},
                "params": {"fasta": fasta, "gtf": gtf},
                "outdir": outdir,
                "resume": False,
            }
            self.assertEqual(
                bulk_salmon_command(context),
                salmon_workflow_command(
                    config,
                    samplesheet=samplesheet,
                    outdir=outdir,
                    fasta=fasta,
                    gtf=gtf,
                ),
            )
            self.assertEqual(
                alevin_command(context),
                alevin_workflow_command(
                    config,
                    samplesheet=samplesheet,
                    outdir=outdir,
                    fasta=fasta,
                    gtf=gtf,
                ),
            )


class SalmonPlanningTests(unittest.TestCase):
    """Native DAG artifacts are known at plan time, unlike the nf-core stages."""

    def _plan(self, preset: str, root: Path):
        target = scaffold_project_preset(preset, root / preset)
        workflow = load_project_config(target / "workflow.toml")
        return plan_workflow(
            workflow,
            load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml")),
        )

    def test_preset_plans_resolve_every_downstream_artifact_path(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            bulk = self._plan("bulk-salmon", root)
            quantify = next(stage for stage in bulk.stages if stage.id == "quantify")
            self.assertEqual(quantify.uses, "bulk.salmon")
            self.assertEqual(
                quantify.outputs["counts"].path,
                quantify.outdir / "counts" / "gene_counts.tsv",
            )
            self.assertEqual(
                quantify.outputs["tx_counts"].path,
                quantify.outdir / "counts" / "transcript_counts.tsv",
            )
            differential = next(
                stage for stage in bulk.stages if stage.id == "differential"
            )
            self.assertIn(
                str(quantify.outputs["counts"].path), " ".join(differential.command)
            )

            single_cell = self._plan("scrnaseq-alevin", root)
            quantify = next(
                stage for stage in single_cell.stages if stage.id == "quantify"
            )
            self.assertEqual(
                quantify.outputs["matrix"].path,
                quantify.outdir / "matrix" / "alevin.h5ad",
            )
            analyze = next(
                stage for stage in single_cell.stages if stage.id == "analyze"
            )
            self.assertIn(
                str(quantify.outputs["matrix"].path), " ".join(analyze.command)
            )


class PackagedSalmonResourceTests(unittest.TestCase):
    def test_native_dags_and_scripts_are_packaged(self) -> None:
        nextflow_root = resources.files("txsuite.resources.nextflow")
        for parts in (("bulk_salmon.nf",), ("single_cell_alevin.nf",), ("modules", "salmon.nf")):
            with self.subTest(name="/".join(parts)):
                self.assertTrue(nextflow_root.joinpath(*parts).is_file())

        salmon_root = resources.files("txsuite.resources.salmon")
        for name in ("Dockerfile", "merge_quants.py"):
            with self.subTest(name=name):
                self.assertTrue(salmon_root.joinpath(name).is_file())

        single_cell_root = resources.files("txsuite.resources.single_cell_python")
        self.assertTrue(single_cell_root.joinpath("alevin_to_h5ad.py").is_file())

    def test_dag_stubs_cover_every_process(self) -> None:
        module = resources.files("txsuite.resources.nextflow").joinpath("modules", "salmon.nf")
        text = module.read_text(encoding="utf-8")
        processes = [line for line in text.splitlines() if line.startswith("process ")]
        self.assertEqual(len(processes), 7)
        self.assertEqual(text.count("stub:"), len(processes))


@unittest.skipIf(
    importlib.util.find_spec("anndata") is None or importlib.util.find_spec("scipy") is None,
    "anndata and scipy are only present inside the single-cell image",
)
class AlevinConversionTests(unittest.TestCase):
    """USA-mode alevin-fry output must become the .h5ad that Scanpy stages read."""

    def setUp(self) -> None:
        self.module = _load_resource_module(
            "txsuite.resources.single_cell_python",
            "alevin_to_h5ad.py",
            "txsuite_alevin_to_h5ad",
        )
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _bundle(self, sample: str, entries: str, *, usa_mode: bool = True) -> Path:
        quant = self.root / sample / "af_quant" / "alevin"
        quant.mkdir(parents=True)
        columns = (
            "gene1-S\ngene2-S\ngene1-U\ngene2-U\ngene1-A\ngene2-A\n"
            if usa_mode
            else "gene1\ngene2\n"
        )
        width = 6 if usa_mode else 2
        rows = entries.strip().splitlines()
        (quant / "quants_mat.mtx").write_text(
            "%%MatrixMarket matrix coordinate real general\n%\n"
            f"2 {width} {len(rows)}\n" + "\n".join(rows) + "\n",
            encoding="utf-8",
        )
        (quant / "quants_mat_rows.txt").write_text(
            "AAACCTGAGAAACCAT\nAAACCTGAGAAACCGC\n", encoding="utf-8"
        )
        (quant / "quants_mat_cols.txt").write_text(columns, encoding="utf-8")
        (quant / "quants_mat.json").write_text(
            '{"usa_mode": %s}\n' % ("true" if usa_mode else "false"), encoding="utf-8"
        )
        return self.root / sample

    def _dense(self, adata):
        matrix = adata.X
        return matrix.toarray() if hasattr(matrix, "toarray") else matrix

    def test_usa_blocks_collapse_to_spliced_plus_ambiguous(self) -> None:
        adata = self.module.load_samples([self._bundle("sample1", "1 1 3\n2 2 5\n1 5 1\n")])

        self.assertEqual(adata.shape, (2, 2))
        self.assertEqual(list(adata.var_names), ["gene1", "gene2"])
        self.assertLessEqual(
            {"counts", "spliced", "unspliced", "ambiguous"},
            {key for key in adata.layers if key},
        )
        dense = self._dense(adata)
        self.assertEqual(dense[0][0], 4)
        self.assertEqual(dense[1][1], 5)
        self.assertEqual(list(adata.obs["sample"].unique()), ["sample1"])

    def test_samples_concatenate_with_unique_barcodes(self) -> None:
        first = self._bundle("sample1", "1 1 3\n")
        second = self._bundle("sample2", "2 2 5\n")

        adata = self.module.load_samples([first, second])

        self.assertEqual(adata.shape[0], 4)
        self.assertEqual(sorted(adata.obs["sample"].unique()), ["sample1", "sample2"])
        self.assertEqual(len(set(adata.obs_names)), 4)

    def test_malformed_bundles_are_rejected(self) -> None:
        with self.assertRaises(self.module.ConversionError):
            self.module.load_samples([self.root / "absent"])
        with self.assertRaises(self.module.ConversionError):
            self.module.load_samples([])

        bundle = self._bundle("sample1", "1 1 3\n")
        (bundle / "af_quant" / "alevin" / "quants_mat_cols.txt").write_text(
            "gene1-S\ngene2-S\ngene1-U\n", encoding="utf-8"
        )
        with self.assertRaises(self.module.ConversionError):
            self.module.load_samples([bundle])


class MergeQuantsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.module = _load_resource_module(
            "txsuite.resources.salmon", "merge_quants.py", "txsuite_merge_quants"
        )
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.tx2gene = self.root / "tx2gene.tsv"
        self.tx2gene.write_text("tx1\tgene1\ntx2\tgene1\ntx3\tgene2\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _quant(self, sample: str, rows: list[tuple[str, float, float]]) -> Path:
        directory = self.root / sample
        directory.mkdir()
        lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        for name, tpm, reads in rows:
            lines.append(f"{name}\t1000\t800.0\t{tpm}\t{reads}")
        (directory / "quant.sf").write_text("\n".join(lines) + "\n", encoding="utf-8")
        return directory

    def test_transcripts_are_summed_per_gene_and_counts_are_integers(self) -> None:
        first = self._quant("sample1", [("tx1", 100.0, 10.4), ("tx2", 50.0, 5.1), ("tx3", 10.0, 1.0)])
        second = self._quant("sample2", [("tx1", 20.0, 2.0), ("tx3", 30.0, 3.0)])

        outputs = self.module.merge([first, second], self.tx2gene, self.root / "counts")

        gene_rows = [
            line.split("\t")
            for line in outputs["gene_counts"].read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(gene_rows[0], ["gene_id", "sample1", "sample2"])
        # 10.4 + 5.1 = 15.5 rounds to 16 under banker-free round-half-away arithmetic
        self.assertEqual(gene_rows[1], ["gene1", "16", "2"])
        self.assertEqual(gene_rows[2], ["gene2", "1", "3"])

        transcript_rows = outputs["transcript_counts"].read_text(encoding="utf-8").splitlines()
        self.assertEqual(transcript_rows[0].split("\t"), ["transcript_id", "sample1", "sample2"])
        self.assertEqual(len(transcript_rows), 4)

    def test_absent_transcripts_are_zero_filled_across_samples(self) -> None:
        first = self._quant("sample1", [("tx1", 100.0, 10.0)])
        second = self._quant("sample2", [("tx3", 30.0, 3.0)])

        outputs = self.module.merge([first, second], self.tx2gene, self.root / "counts")
        rows = dict(
            (line.split("\t")[0], line.split("\t")[1:])
            for line in outputs["gene_counts"].read_text(encoding="utf-8").splitlines()[1:]
        )
        self.assertEqual(rows["gene1"], ["10", "0"])
        self.assertEqual(rows["gene2"], ["0", "3"])

    def test_inconsistent_or_unmatched_inputs_fail_loudly(self) -> None:
        quant = self._quant("sample1", [("unknown_tx", 1.0, 1.0)])
        with self.assertRaises(self.module.MergeError):
            self.module.merge([quant], self.tx2gene, self.root / "counts")

        conflicting = self.root / "conflicting.tsv"
        conflicting.write_text("tx1\tgene1\ntx1\tgene2\n", encoding="utf-8")
        with self.assertRaises(self.module.MergeError):
            self.module.merge([quant], conflicting, self.root / "counts")

        with self.assertRaises(self.module.MergeError):
            self.module.merge([self.root / "absent"], self.tx2gene, self.root / "counts")


if __name__ == "__main__":
    unittest.main()
