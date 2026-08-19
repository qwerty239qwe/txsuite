"""Native STAR genome alignment and gene counting.

STAR emits gene counts alongside the alignment, so the delicate part is not the
alignment call but choosing the right count column: ReadsPerGene.out.tab holds
unstranded, forward, and reverse totals, and picking wrong produces a matrix
that is mostly noise without failing. Those tests run the real inference.
"""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.alignment import STRANDEDNESS, align_workflow_command
from txsuite.config import DEFAULT_CONFIG, load_config
from txsuite.project.adapters.bulk import bulk_align_command
from txsuite.project.config import load_project_config
from txsuite.project.planner import PlanningError, plan_workflow
from txsuite.project.presets import scaffold_project_preset
from txsuite.project.registry import get_stage_spec, validate_stage_parameters
from txsuite.runtime import TxSuiteError


def _merge_module():
    path = Path(
        str(resources.files("txsuite.resources.star").joinpath("merge_star_counts.py"))
    )
    spec = importlib.util.spec_from_file_location("txsuite_merge_star_counts", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AlignCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.samplesheet = self.root / "samplesheet.csv"
        self.samplesheet.write_text(
            "sample,fastq_1,fastq_2\nsample1,r1.fastq.gz,r2.fastq.gz\n", encoding="utf-8"
        )
        self.index = self.root / "star_index"
        self.index.mkdir()
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _command(self, **kwargs):
        arguments = {
            "samplesheet": self.samplesheet,
            "star_index": self.index,
            "outdir": self.outdir,
        }
        arguments.update(kwargs)
        return align_workflow_command(self.config, **arguments)

    def test_command_names_the_dag_index_and_owned_image(self) -> None:
        command = self._command()

        self.assertTrue(command[2].endswith("bulk_align.nf"))
        self.assertEqual(
            command[command.index("--star_index") + 1], str(self.index.resolve())
        )
        self.assertEqual(
            command[command.index("--star_image") + 1], self.config["images"]["star"]
        )
        self.assertEqual(command[command.index("--star_strandedness") + 1], "auto")
        self.assertEqual(command[command.index("--star_two_pass") + 1], "true")
        self.assertNotIn("-resume", command)

    def test_strandedness_and_two_pass_are_forwarded(self) -> None:
        for value in STRANDEDNESS:
            with self.subTest(strandedness=value):
                command = self._command(strandedness=value)
                self.assertEqual(
                    command[command.index("--star_strandedness") + 1], value
                )
        self.assertEqual(
            self._command(two_pass=False)[
                self._command(two_pass=False).index("--star_two_pass") + 1
            ],
            "false",
        )

    def test_invalid_settings_are_rejected(self) -> None:
        for kwargs in (
            {"strandedness": "sideways"},
            {"threads": 0},
            {"memory_gb": 0},
            {"star_index": self.root / "absent"},
            {"samplesheet": self.root / "absent.csv"},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                self._command(**kwargs)

    def test_adapter_matches_the_library_helper(self) -> None:
        context = {
            "global_config": self.config,
            "resolved_inputs": {
                "samplesheet": self.samplesheet,
                "star_index": self.index,
            },
            "params": {"strandedness": "reverse"},
            "outdir": self.outdir,
        }
        self.assertEqual(
            bulk_align_command(context), self._command(strandedness="reverse")
        )


class AlignStageTests(unittest.TestCase):
    def test_stage_produces_counts_that_feed_differential_expression(self) -> None:
        spec = get_stage_spec("bulk.align")

        self.assertEqual(spec.modality, "bulk")
        self.assertEqual(set(spec.inputs), {"samplesheet", "star_index"})
        self.assertEqual(
            spec.inputs["star_index"],
            get_stage_spec("bulk.star-reference").outputs["star_index"],
        )
        self.assertEqual(
            spec.outputs["counts"], get_stage_spec("bulk.de").inputs["counts"]
        )
        self.assertEqual(spec.output_policies["alignments"].kind, "directory")
        self.assertEqual(spec.required_images, ("images.star",))
        self.assertTrue(spec.supports_resume)

    def test_samplesheet_type_is_weaker_than_the_nfcore_contract(self) -> None:
        # bulk.rnaseq requires a strandedness column; the native stages infer it,
        # so wiring a native sheet into bulk.rnaseq must not type-check.
        self.assertEqual(
            get_stage_spec("bulk.align").inputs["samplesheet"],
            "bulk.fastq-samplesheet",
        )
        self.assertEqual(
            get_stage_spec("bulk.salmon").inputs["samplesheet"],
            "bulk.fastq-samplesheet",
        )
        self.assertNotEqual(
            get_stage_spec("bulk.rnaseq").inputs["samplesheet"],
            "bulk.fastq-samplesheet",
        )

    def test_defaults_and_validation(self) -> None:
        spec = get_stage_spec("bulk.align")
        defaults = validate_stage_parameters(spec, {})
        self.assertEqual(defaults["strandedness"], "auto")
        self.assertTrue(defaults["two_pass"])
        for params in ({"strandedness": "sideways"}, {"threads": 0}):
            with self.subTest(params=params), self.assertRaises(TxSuiteError):
                validate_stage_parameters(spec, params)


class AlignPlanningTests(unittest.TestCase):
    def _plan(self):
        target = Path(tempfile.mkdtemp()) / "project"
        scaffold_project_preset("bulk-align", target)
        return plan_workflow(load_project_config(target / "workflow.toml"), DEFAULT_CONFIG)

    def test_preset_wires_index_and_counts_end_to_end(self) -> None:
        plan = self._plan()
        reference, align, de = (plan.stage(i) for i in ("ref", "align", "differential"))

        self.assertIn("ref", align.depends_on)
        self.assertEqual(
            align.inputs["star_index"].value, reference.outputs["star_index"].path
        )
        self.assertEqual(
            align.outputs["counts"].path, align.outdir / "counts" / "gene_counts.tsv"
        )
        self.assertEqual(
            align.outputs["strandedness"].path,
            align.outdir / "counts" / "strandedness.tsv",
        )
        # The whole point of emitting counts: differential expression consumes
        # them without an intervening quantification stage.
        self.assertIn("align", de.depends_on)
        self.assertIn(str(align.outputs["counts"].path), " ".join(de.command))

    def test_alignment_artifacts_resolve_at_plan_time(self) -> None:
        align = self._plan().stage("align")
        for name in ("alignments", "logs", "counts", "strandedness", "results"):
            with self.subTest(output=name):
                self.assertIsNotNone(align.outputs[name].path)
        self.assertEqual(dict(align.postflight), {})

    def test_index_is_required(self) -> None:
        from txsuite.project.config import parse_project_config

        workflow = parse_project_config(
            {
                "schema_version": 1,
                "project": {"id": "a", "modality": "bulk", "output_root": "results"},
                "execution": {"profile": "docker", "resume": True},
                "workflow": {
                    "stages": [
                        {
                            "id": "align",
                            "uses": "bulk.align",
                            "inputs": {"samplesheet": "samplesheet.csv"},
                        }
                    ]
                },
            },
            source_path=Path("/work/workflow.toml"),
        )
        with self.assertRaisesRegex(PlanningError, "missing required input"):
            plan_workflow(workflow, DEFAULT_CONFIG)


class StrandednessInferenceTests(unittest.TestCase):
    """The unstranded column is roughly forward + reverse, not a competitor.

    A naive "largest column wins" rule therefore always answers unstranded. The
    signal is the split between the two stranded columns.
    """

    def setUp(self) -> None:
        self.module = _merge_module()
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _sample(self, name: str, rows: dict[str, tuple[int, int, int]]) -> Path:
        directory = self.root / name
        directory.mkdir(parents=True, exist_ok=True)
        lines = [
            "N_unmapped\t100\t100\t100",
            "N_multimapping\t50\t50\t50",
            # Large enough to invert the answer if summary rows leaked into the
            # totals, which is exactly the mistake worth guarding.
            "N_noFeature\t900000\t900000\t900000",
            "N_ambiguous\t10\t10\t10",
        ]
        lines += ["\t".join([gene, *map(str, values)]) for gene, values in rows.items()]
        (directory / "ReadsPerGene.out.tab").write_text(
            "\n".join(lines) + "\n", encoding="utf-8"
        )
        return directory

    def _counts(self, outdir: Path) -> list[list[str]]:
        return [
            line.split("\t")
            for line in (outdir / "gene_counts.tsv").read_text(encoding="utf-8").splitlines()
        ]

    def test_each_orientation_selects_its_own_column(self) -> None:
        cases = {
            "reverse": ({"g1": (1000, 20, 980)}, "980"),
            "forward": ({"g1": (1000, 980, 20)}, "980"),
            "unstranded": ({"g1": (1000, 505, 495)}, "1000"),
        }
        for expected, (rows, value) in cases.items():
            with self.subTest(strandedness=expected):
                outdir = self.root / f"out_{expected}"
                self.module.merge([self._sample(expected, rows)], outdir)
                report = (outdir / "strandedness.tsv").read_text(
                    encoding="utf-8"
                ).splitlines()[1].split("\t")
                self.assertEqual(report[1], expected)
                self.assertEqual(report[2], "inferred")
                self.assertEqual(self._counts(outdir)[1][1], value)

    def test_summary_rows_are_excluded_from_the_matrix_and_the_totals(self) -> None:
        outdir = self.root / "out_summary"
        self.module.merge([self._sample("s", {"g1": (1000, 20, 980)})], outdir)
        rows = self._counts(outdir)
        self.assertEqual([row[0] for row in rows[1:]], ["g1"])
        report = (outdir / "strandedness.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        self.assertEqual(report[1], "reverse")
        self.assertEqual(report[4], "20")
        self.assertEqual(report[5], "980")

    def test_an_unclear_split_is_refused_rather_than_guessed(self) -> None:
        with self.assertRaisesRegex(self.module.MergeError, "cannot infer strandedness"):
            self.module.merge(
                [self._sample("ambiguous", {"g1": (1000, 700, 300)})],
                self.root / "out_ambiguous",
            )

    def test_samples_that_disagree_are_refused(self) -> None:
        with self.assertRaisesRegex(self.module.MergeError, "disagree on strandedness"):
            self.module.merge(
                [
                    self._sample("a", {"g1": (1000, 20, 980)}),
                    self._sample("b", {"g1": (1000, 980, 20)}),
                ],
                self.root / "out_disagree",
            )

    def test_no_assigned_reads_is_refused(self) -> None:
        with self.assertRaises(self.module.MergeError):
            self.module.merge(
                [self._sample("empty", {"g1": (0, 0, 0)})], self.root / "out_empty"
            )

    def test_explicit_strandedness_overrides_inference_and_is_recorded(self) -> None:
        outdir = self.root / "out_override"
        self.module.merge(
            [self._sample("s", {"g1": (1000, 20, 980)})], outdir, strandedness="forward"
        )
        self.assertEqual(self._counts(outdir)[1][1], "20")
        report = (outdir / "strandedness.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        self.assertEqual(report[1], "forward")
        self.assertEqual(report[2], "explicit")

    def test_genes_absent_from_one_sample_are_zero_filled(self) -> None:
        outdir = self.root / "out_union"
        self.module.merge(
            [
                self._sample("a", {"g1": (1000, 20, 980)}),
                self._sample("b", {"g2": (500, 10, 490)}),
            ],
            outdir,
        )
        rows = {row[0]: row[1:] for row in self._counts(outdir)[1:]}
        self.assertEqual(rows["g1"], ["980", "0"])
        self.assertEqual(rows["g2"], ["0", "490"])

    def test_malformed_inputs_are_refused(self) -> None:
        broken = self.root / "broken"
        broken.mkdir()
        (broken / "ReadsPerGene.out.tab").write_text("g1\t1\n", encoding="utf-8")
        with self.assertRaises(self.module.MergeError):
            self.module.merge([broken], self.root / "out_broken")
        with self.assertRaises(self.module.MergeError):
            self.module.merge([self.root / "absent"], self.root / "out_absent")
        with self.assertRaises(self.module.MergeError):
            self.module.merge([], self.root / "out_none")


class PackagedAlignResourceTests(unittest.TestCase):
    def test_dag_and_merge_script_are_packaged(self) -> None:
        nextflow_root = resources.files("txsuite.resources.nextflow")
        self.assertTrue(nextflow_root.joinpath("bulk_align.nf").is_file())
        module = nextflow_root.joinpath("modules", "star_align.nf")
        self.assertTrue(module.is_file())
        text = module.read_text(encoding="utf-8")
        processes = [line for line in text.splitlines() if line.startswith("process ")]
        self.assertEqual(len(processes), 3)
        self.assertEqual(text.count("stub:"), len(processes))
        self.assertIn("--quantMode GeneCounts", text)
        # Each sample's BAM has the same STAR-assigned filename, so both the BAM
        # and its index must publish under a per-sample directory.
        self.assertIn(
            'publishDir "${params.outdir}/alignments/${sample}", mode: \'copy\'', text
        )
        self.assertIn("pattern: '*/*.bam'", text)
        self.assertIn("pattern: '*/*.out'", text)
        self.assertTrue(
            resources.files("txsuite.resources.star")
            .joinpath("merge_star_counts.py")
            .is_file()
        )

    def test_image_recipe_ships_the_merge_script(self) -> None:
        dockerfile = (
            resources.files("txsuite.resources.star")
            .joinpath("Dockerfile")
            .read_text(encoding="utf-8")
        )
        self.assertIn("COPY merge_star_counts.py /opt/txsuite/merge_star_counts.py", dockerfile)


if __name__ == "__main__":
    unittest.main()
