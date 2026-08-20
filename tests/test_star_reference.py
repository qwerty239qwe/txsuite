"""Shared STAR/GATK reference preparation.

A STAR index costs about an hour and 32 GB for a human genome, so it is built
once as its own stage and reused. These tests cover the command contract, the
overhang derivation that silently degrades results when wrong, and the packaged
DAG resources.
"""

from __future__ import annotations

import tempfile
import unittest
from contextlib import redirect_stdout
from importlib import resources
from io import StringIO
from pathlib import Path

from txsuite.cli import run
from txsuite.config import DEFAULT_CONFIG, load_config
from txsuite.project.adapters.bulk import bulk_star_reference_command
from txsuite.project.config import parse_project_config
from txsuite.project.planner import PlanningError, plan_workflow
from txsuite.project.registry import (
    get_stage_spec,
    list_stage_specs,
    validate_stage_parameters,
)
from txsuite.reference import build_star_image, star_reference_workflow_command
from txsuite.runtime import TxSuiteError


class StarReferenceCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.fasta = self.root / "genome.fa"
        self.fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.gtf = self.root / "genes.gtf"
        self.gtf.write_text("#gtf\n", encoding="utf-8")
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _command(self, **kwargs):
        arguments = {"fasta": self.fasta, "gtf": self.gtf, "outdir": self.outdir}
        arguments.update(kwargs)
        return star_reference_workflow_command(self.config, **arguments)

    def test_overhang_defaults_to_read_length_minus_one(self) -> None:
        command = self._command(read_length=150)
        self.assertEqual(command[command.index("--star_read_length") + 1], "150")
        self.assertEqual(command[command.index("--star_sjdb_overhang") + 1], "149")

    def test_explicit_overhang_overrides_the_derived_value(self) -> None:
        command = self._command(read_length=150, sjdb_overhang=100)
        self.assertEqual(command[command.index("--star_sjdb_overhang") + 1], "100")

    def test_command_names_the_owned_image_and_native_dag(self) -> None:
        command = self._command()
        self.assertTrue(command[2].endswith("bulk_star_reference.nf"))
        self.assertEqual(
            command[command.index("--star_image") + 1], self.config["images"]["star"]
        )
        self.assertNotIn("--star_sa_index_nbases", command)
        self.assertNotIn("-resume", command)
        self.assertIn("-resume", self._command(resume=True))

    def test_invalid_sizes_are_rejected(self) -> None:
        for kwargs in (
            {"read_length": 1},
            {"sjdb_overhang": 0},
            {"threads": 0},
            {"memory_gb": 0},
            {"genome_sa_index_nbases": 0},
            {"genome_sa_index_nbases": 17},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                self._command(**kwargs)

    def test_missing_reference_files_are_rejected(self) -> None:
        with self.assertRaises(TxSuiteError):
            self._command(fasta=self.root / "absent.fa")
        with self.assertRaises(TxSuiteError):
            self._command(gtf=self.root / "absent.gtf")

    def test_adapter_matches_the_library_helper(self) -> None:
        context = {
            "global_config": self.config,
            "resolved_inputs": {"fasta": self.fasta, "gtf": self.gtf},
            "params": {"read_length": 75},
            "outdir": self.outdir,
        }
        self.assertEqual(
            bulk_star_reference_command(context),
            self._command(read_length=75),
        )


class StarReferenceStageTests(unittest.TestCase):
    def test_stage_declares_every_reusable_artifact(self) -> None:
        spec = get_stage_spec("bulk.star-reference")

        self.assertEqual(set(spec.inputs), {"fasta", "gtf"})
        self.assertEqual(
            set(spec.outputs),
            {"results", "star_index", "fasta_fai", "dict", "manifest"},
        )
        self.assertEqual(spec.output_policies["star_index"].kind, "directory")
        self.assertEqual(spec.required_images, ("images.star",))
        self.assertEqual(dict(spec.optional_inputs), {})

    def test_defaults_and_validation(self) -> None:
        spec = get_stage_spec("bulk.star-reference")
        defaults = validate_stage_parameters(spec, {})
        self.assertEqual(defaults["read_length"], 100)
        self.assertIsNone(defaults["sjdb_overhang"])
        for params in ({"read_length": 0}, {"threads": -1}, {"memory_gb": 0}):
            with self.subTest(params=params), self.assertRaises(TxSuiteError):
                validate_stage_parameters(spec, params)


class ReferenceBasenameTests(unittest.TestCase):
    """GATK resolves the index and dictionary from the reference basename.

    Publishing them under a canonical name made the planner's job easy and made
    the artifacts unusable: GATK given `--reference GRCh38.fa` looks for
    `GRCh38.fa.fai` and `GRCh38.dict`, not `genome.*`. The preset hid this by
    happening to name its placeholder `genome.fa`.
    """

    def _plan(self, fasta_name: str):
        workflow = parse_project_config(
            {
                "schema_version": 1,
                "project": {
                    "id": "reference",
                    "modality": "bulk",
                    "output_root": "results",
                },
                "execution": {"profile": "docker", "resume": True},
                "workflow": {
                    "stages": [
                        {
                            "id": "ref",
                            "uses": "bulk.star-reference",
                            "inputs": {"fasta": fasta_name, "gtf": "genes.gtf"},
                        }
                    ]
                },
            },
            source_path=Path("/work/workflow.toml"),
        )
        return plan_workflow(workflow, DEFAULT_CONFIG).stage("ref")

    def test_index_names_follow_the_callers_fasta(self) -> None:
        for fasta_name, fai, dictionary in (
            ("genome.fa", "genome.fa.fai", "genome.dict"),
            ("GRCh38.primary.fa", "GRCh38.primary.fa.fai", "GRCh38.primary.dict"),
            ("data/hg38.fasta", "hg38.fasta.fai", "hg38.dict"),
        ):
            with self.subTest(fasta=fasta_name):
                reference = self._plan(fasta_name)
                self.assertEqual(
                    reference.outputs["fasta_fai"].path,
                    reference.outdir / "reference" / fai,
                )
                self.assertEqual(
                    reference.outputs["dict"].path,
                    reference.outdir / "reference" / dictionary,
                )

    def test_dict_name_drops_only_the_final_extension(self) -> None:
        # GATK wants X.dict for X.fa; a multi-dot name must keep its inner dots.
        reference = self._plan("Homo_sapiens.GRCh38.dna.fa")
        self.assertEqual(
            reference.outputs["dict"].path.name, "Homo_sapiens.GRCh38.dna.dict"
        )


class PackagedStarResourceTests(unittest.TestCase):
    def test_dag_and_image_recipe_are_packaged(self) -> None:
        nextflow_root = resources.files("txsuite.resources.nextflow")
        self.assertTrue(nextflow_root.joinpath("bulk_star_reference.nf").is_file())
        module = nextflow_root.joinpath("modules", "star_reference.nf")
        self.assertTrue(module.is_file())

        text = module.read_text(encoding="utf-8")
        processes = [line for line in text.splitlines() if line.startswith("process ")]
        self.assertEqual(len(processes), 4)
        self.assertEqual(text.count("stub:"), len(processes))
        # GATK resolves the index and dictionary from the reference basename, so
        # the DAG must keep the caller's FASTA name rather than a canonical one.
        self.assertNotIn("stageAs: 'genome.fa'", text)
        self.assertIn('path "${fasta}.fai", emit: fai', text)
        self.assertIn('path "${fasta.baseName}.dict", emit: dict', text)

        self.assertTrue(
            resources.files("txsuite.resources.star").joinpath("Dockerfile").is_file()
        )

    def test_dag_avoids_leading_operator_continuations(self) -> None:
        # Groovy terminates a statement at a line break when the expression is
        # already complete, so a ternary split across lines is a syntax error
        # unless the previous line ends with an explicit backslash.
        root = resources.files("txsuite.resources.nextflow")
        scripts = [child for child in root.iterdir() if child.name.endswith(".nf")]
        scripts += [
            child
            for child in root.joinpath("modules").iterdir()
            if child.name.endswith(".nf")
        ]
        self.assertGreaterEqual(len(scripts), 8)
        for script in scripts:
            with self.subTest(name=script.name):
                lines = script.read_text(encoding="utf-8").splitlines()
                for number, line in enumerate(lines, start=1):
                    stripped = line.strip()
                    if not (stripped.startswith("?") or stripped.startswith(":")):
                        continue
                    previous = lines[number - 2].rstrip() if number > 1 else ""
                    self.assertTrue(
                        previous.endswith("\\"),
                        f"{script.name} line {number} continues a ternary without an "
                        f"explicit backslash: {stripped}",
                    )



class ExecutionProfileTests(unittest.TestCase):
    """Which stages may run under Apptainer.

    Profile support used to be a hardcoded list of three stage IDs, so every
    Nextflow DAG added afterwards was silently docker-only even though Nextflow
    picks the runtime through its own profile. It is now derived from the
    executables a stage declares.
    """

    def _plan(self, uses: str, profile: str):
        stages = {
            "bulk.star-reference": {
                "inputs": {"fasta": "genome.fa", "gtf": "genes.gtf"},
                "params": {},
            },
            "bulk.salmon": {
                "inputs": {"samplesheet": "samplesheet.csv"},
                "params": {"fasta": "genome.fa", "gtf": "genes.gtf"},
            },
            "bulk.de": {
                "inputs": {"counts": "counts.tsv", "metadata": "metadata.tsv"},
                "params": {
                    "design": "condition",
                    "reference": "control",
                    "test": "treated",
                },
            },
        }[uses]
        workflow = parse_project_config(
            {
                "schema_version": 1,
                "project": {
                    "id": "profiles",
                    "modality": "bulk",
                    "output_root": "results",
                },
                "execution": {"profile": profile, "resume": True},
                "workflow": {
                    "stages": [{"id": "stage", "uses": uses, **stages}]
                },
            },
            source_path=Path("/work/workflow.toml"),
        )
        return plan_workflow(workflow, DEFAULT_CONFIG)

    def test_nextflow_stages_support_both_profiles(self) -> None:
        for uses in ("bulk.star-reference", "bulk.salmon"):
            for profile in ("docker", "apptainer"):
                with self.subTest(uses=uses, profile=profile):
                    plan = self._plan(uses, profile)
                    command = plan.stage("stage").command
                    self.assertEqual(command[command.index("-profile") + 1], profile)

    def test_docker_only_stages_still_reject_apptainer(self) -> None:
        self._plan("bulk.de", "docker")
        with self.assertRaises(PlanningError):
            self._plan("bulk.de", "apptainer")

    def test_every_nextflow_stage_declares_the_executable(self) -> None:
        # The profile rule reads required_executables, so a Nextflow stage that
        # forgets to declare it would quietly lose Apptainer support.
        nextflow_stages = {
            "bulk.rnaseq",
            "bulk.rnavar",
            "bulk.salmon",
            "bulk.star-reference",
            "single-cell.scrnaseq",
            "single-cell.alevin",
            "single-cell.pseudobulk-de",
        }
        covered = {spec.uses for spec in list_stage_specs()} & nextflow_stages
        self.assertEqual(covered, nextflow_stages)
        for spec in list_stage_specs():
            if spec.uses not in nextflow_stages:
                continue
            with self.subTest(uses=spec.uses):
                self.assertIn("nextflow", spec.required_executables)


class StarImageBuildTests(unittest.TestCase):
    def test_empty_and_malformed_tags_are_rejected_before_docker_runs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory) / "run"
            for tag in ("", "   ", "bad tag", "-leading-dash"):
                with self.subTest(tag=tag), self.assertRaises(TxSuiteError):
                    build_star_image(tag, run_dir=run_dir)
            self.assertFalse(run_dir.exists())

    def test_env_build_star_is_wired_into_the_cli(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            status = run(
                ["env", "build", "star", "--config", "txsuite.toml", "--dry-run"]
            )
        self.assertEqual(status, 0)
        printed = output.getvalue()
        self.assertIn("docker build --tag txsuite/star:0.2.0", printed)
        self.assertIn("<bundled-star-context>", printed)

    def test_recipe_pins_every_tool_the_reference_dag_calls(self) -> None:
        dockerfile = (
            resources.files("txsuite.resources.star")
            .joinpath("Dockerfile")
            .read_text(encoding="utf-8")
        )
        for pin in ("star=2.7.11b", "samtools=1.24", "gatk4=4.6.2.0"):
            with self.subTest(pin=pin):
                self.assertIn(pin, dockerfile)
        self.assertIn("@sha256:", dockerfile.splitlines()[0])


if __name__ == "__main__":
    unittest.main()
