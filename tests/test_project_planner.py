from __future__ import annotations

import copy
import json
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path
from unittest.mock import patch

from txsuite.config import DEFAULT_CONFIG
from txsuite.project.config import parse_project_config
from txsuite.project.model import (
    ArtifactReference,
    ExecutionSettings,
    ProjectDefinition,
    ResolvedWorkflow,
    StageDefinition,
)
from txsuite.project.plan_format import format_plan_human, render_plan
from txsuite.project.planner import PlanningError, plan_workflow
from txsuite.project.results import verify_required_outputs


SOURCE = Path("/work/project.toml")
FIXTURES = Path(__file__).parent / "fixtures" / "plans"


def resolved(
    modality: str,
    stages: list[dict[str, object]],
    *,
    profile: str = "docker",
) -> ResolvedWorkflow:
    return parse_project_config(
        {
            "schema_version": 1,
            "project": {
                "id": "example",
                "modality": modality,
                "output_root": "results",
            },
            "execution": {"profile": profile, "resume": True},
            "workflow": {"stages": stages},
        },
        source_path=SOURCE,
    )


def summary(plan) -> dict[str, object]:
    return {
        "order": [stage.id for stage in plan.stages],
        "states": {stage.id: stage.command_state for stage in plan.stages},
        "dependencies": {
            stage.id: list(stage.depends_on) for stage in plan.stages
        },
        "outputs": {
            stage.id: {
                name: (output.value.as_posix().removeprefix(output.value.drive)
                       if isinstance(output.value, Path) else str(output.value))
                for name, output in stage.outputs.items()
            }
            for stage in plan.stages
        },
    }


class ProjectPlannerGoldenTests(unittest.TestCase):
    def assert_golden(self, name: str, plan) -> None:
        expected = json.loads((FIXTURES / f"{name}.json").read_text(encoding="utf-8"))
        self.assertEqual(summary(plan), expected)

    def test_bulk_golden_resolves_owned_outputs_and_defers_nfcore_artifacts(self) -> None:
        workflow = resolved(
            "bulk",
            [
                {
                    "id": "raw",
                    "uses": "bulk.rnaseq",
                    "inputs": {"samplesheet": "samples.csv"},
                },
                {
                    "id": "de",
                    "uses": "bulk.de",
                    "inputs": {
                        "counts": "${raw.counts}",
                        "metadata": "metadata.tsv",
                    },
                    "params": {
                        "design": "condition",
                        "reference": "control",
                        "test": "treated",
                    },
                },
                {
                    "id": "enrich",
                    "uses": "bulk.enrichment",
                    "inputs": {
                        "de_results": "${de.de_results}",
                        "genesets": "sets.gmt",
                    },
                },
            ],
        )
        plan = plan_workflow(workflow, DEFAULT_CONFIG)
        self.assert_golden("bulk", plan)
        self.assertIn("${raw.counts}", " ".join(plan.stage("de").command))
        self.assertEqual(plan.stage("de").params["method"], "deseq2")

    def test_scrna_golden(self) -> None:
        plan = plan_workflow(
            resolved(
                "single-cell",
                [
                    {
                        "id": "raw",
                        "uses": "single-cell.scrnaseq",
                        "inputs": {"samplesheet": "samples.csv"},
                    },
                    {
                        "id": "scanpy",
                        "uses": "single-cell.scanpy",
                        "inputs": {"input": "${raw.matrix}"},
                    },
                ],
            ),
            DEFAULT_CONFIG,
        )
        self.assert_golden("scrna", plan)
        raw = plan.stage("raw")
        self.assertEqual(raw.outputs["matrix"].kind, "file")
        self.assertTrue(raw.outputs["matrix"].non_empty)
        self.assertEqual(
            dict(raw.postflight["policies"]["matrix"]),
            {"kind": "file", "non_empty": True},
        )
        scanpy = plan.stage("scanpy")
        self.assertIn(
            "type=bind,source=${raw.matrix},target=/input/data.h5ad,readonly",
            scanpy.command,
        )
        self.assertIn("/input/data.h5ad", scanpy.command)

    def test_pseudobulk_golden(self) -> None:
        plan = plan_workflow(
            resolved(
                "single-cell",
                [
                    {
                        "id": "pseudobulk",
                        "uses": "single-cell.pseudobulk",
                        "inputs": {"h5ad": "analysis.h5ad"},
                        "params": {"sample_column": "sample", "design": "condition"},
                    },
                    {
                        "id": "native_de",
                        "uses": "single-cell.pseudobulk-de",
                        "inputs": {"h5ad": "analysis.h5ad"},
                        "params": {
                            "sample_column": "sample",
                            "design": "condition",
                            "reference": "control",
                            "test": "treated",
                        },
                    },
                ],
            ),
            DEFAULT_CONFIG,
        )
        self.assert_golden("pseudobulk", plan)

    def test_branching_golden_is_independent_of_declaration_order(self) -> None:
        stages = [
            {
                "id": "z_de",
                "uses": "bulk.de",
                "inputs": {"counts": "${raw.counts}", "metadata": "z.tsv"},
                "params": {"design": "group", "reference": "a", "test": "b"},
            },
            {
                "id": "a_de",
                "uses": "bulk.de",
                "inputs": {"counts": "${raw.counts}", "metadata": "a.tsv"},
                "params": {"design": "group", "reference": "a", "test": "b"},
            },
            {
                "id": "raw",
                "uses": "bulk.rnaseq",
                "inputs": {"samplesheet": "samples.csv"},
                "outputs": {"counts": "artifacts/counts.tsv"},
            },
        ]
        plan = plan_workflow(resolved("bulk", stages), DEFAULT_CONFIG)
        self.assert_golden("branching", plan)
        reordered = plan_workflow(resolved("bulk", list(reversed(stages))), DEFAULT_CONFIG)
        self.assertEqual(plan.to_dict(), reordered.to_dict())


class ProjectPlannerValidationTests(unittest.TestCase):
    def test_relative_path_params_are_anchored_to_workflow_source_not_cwd(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            project_dir = root / "project"
            config_dir = project_dir / "config"
            config_dir.mkdir(parents=True)
            params_file = config_dir / "params.json"
            params_file.write_text("{}\n", encoding="utf-8")
            workflow = parse_project_config(
                {
                    "schema_version": 1,
                    "project": {
                        "id": "cwd_independent",
                        "modality": "bulk",
                        "output_root": "results",
                    },
                    "execution": {"profile": "docker", "resume": False},
                    "workflow": {
                        "stages": [{
                            "id": "raw",
                            "uses": "bulk.rnaseq",
                            "inputs": {"samplesheet": "samples.csv"},
                            "params": {"params_file": "config/params.json"},
                        }]
                    },
                },
                source_path=project_dir / "workflow.toml",
            )

            with patch("os.getcwd", return_value=str(root / "cwd-one")):
                first = plan_workflow(workflow, DEFAULT_CONFIG)
            with patch("os.getcwd", return_value=str(root / "cwd-two")):
                second = plan_workflow(workflow, DEFAULT_CONFIG)

            expected = params_file.resolve()
            self.assertEqual(first.stage("raw").params["params_file"], expected)
            self.assertEqual(first.stage("raw").command, second.stage("raw").command)
            option = first.stage("raw").command.index("-params-file")
            self.assertEqual(first.stage("raw").command[option + 1], str(expected))

    def test_project_profile_is_authoritative_and_backend_compatible(self) -> None:
        global_config = copy.deepcopy(DEFAULT_CONFIG)
        global_config["execution"]["profile"] = "apptainer"
        workflow = resolved(
            "bulk",
            [{"id": "raw", "uses": "bulk.rnaseq", "inputs": {"samplesheet": "s.csv"}}],
            profile="docker",
        )
        plan = plan_workflow(workflow, global_config)
        command = plan.stage("raw").command
        self.assertEqual(command[command.index("-profile") + 1], "docker")
        self.assertEqual(global_config["execution"]["profile"], "apptainer")

        raw_apptainer = plan_workflow(
            resolved(
                "bulk",
                [{"id": "raw", "uses": "bulk.rnaseq", "inputs": {"samplesheet": "s.csv"}}],
                profile="apptainer",
            ),
            DEFAULT_CONFIG,
        )
        self.assertIn("apptainer", raw_apptainer.stage("raw").command)

        native = plan_workflow(
            resolved(
                "single-cell",
                [{
                    "id": "native",
                    "uses": "single-cell.pseudobulk-de",
                    "inputs": {"h5ad": "a.h5ad"},
                    "params": {
                        "sample_column": "sample",
                        "design": "condition",
                        "reference": "a",
                        "test": "b",
                    },
                }],
                profile="apptainer",
            ),
            DEFAULT_CONFIG,
        )
        self.assertIn("apptainer", native.stage("native").command)

        direct = resolved(
            "bulk",
            [{
                "id": "de",
                "uses": "bulk.de",
                "inputs": {"counts": "c.tsv", "metadata": "m.tsv"},
                "params": {"design": "group", "reference": "a", "test": "b"},
            }],
            profile="apptainer",
        )
        with self.assertRaisesRegex(PlanningError, "currently emits Docker commands"):
            plan_workflow(direct, DEFAULT_CONFIG)

    def test_raw_nfcore_postflight_metadata_is_canonical(self) -> None:
        plan = plan_workflow(
            resolved(
                "bulk",
                [{
                    "id": "raw",
                    "uses": "bulk.rnaseq",
                    "inputs": {"samplesheet": "s.csv"},
                    "outputs": {"counts": "artifacts/counts.tsv"},
                }],
            ),
            DEFAULT_CONFIG,
        )
        postflight = plan.stage("raw").postflight
        self.assertEqual(postflight["adapter"], "nfcore")
        self.assertEqual(postflight["pipeline"], "nf-core/rnaseq")
        self.assertEqual(postflight["release"], "3.26.0")
        self.assertEqual(postflight["results_root"], Path("/work/results/raw").resolve())
        self.assertEqual(postflight["artifacts"]["counts"], "bulk.gene-counts")
        self.assertEqual(
            dict(postflight["policies"]["counts"]),
            {"kind": "file", "non_empty": True},
        )
        self.assertEqual(postflight["overrides"]["counts"], Path("/work/artifacts/counts.tsv").resolve())

    def test_planned_output_policies_reject_wrong_kind_zero_byte_and_empty_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scanpy_workflow = parse_project_config(
                {
                    "schema_version": 1,
                    "project": {"id": "verify", "modality": "single-cell", "output_root": "results"},
                    "execution": {"profile": "docker", "resume": False},
                    "workflow": {"stages": [{"id": "scan", "uses": "single-cell.scanpy", "inputs": {"input": "matrix"}}]},
                },
                source_path=root / "workflow.toml",
            )
            scanpy = plan_workflow(scanpy_workflow, DEFAULT_CONFIG).stage("scan")
            required = scanpy.required_outputs
            expected_output = (root / "results" / "scan" / "analysis.h5ad").resolve()
            self.assertEqual(len(required), 1)
            self.assertEqual(
                dict(required[0]),
                {
                    "path": expected_output,
                    "kind": "file",
                    "non_empty": True,
                },
            )
            self.assertEqual(
                scanpy.to_dict()["required_outputs"],
                [{
                    "path": str(expected_output),
                    "kind": "file",
                    "non_empty": True,
                }],
            )

            output = required[0]["path"]
            output.parent.mkdir(parents=True)
            output.mkdir()
            self.assertFalse(verify_required_outputs(required).ok)
            output.rmdir()
            output.write_bytes(b"")
            self.assertFalse(verify_required_outputs(required).ok)
            output.write_bytes(b"h5ad")
            self.assertTrue(verify_required_outputs(required).ok)

            enrichment_workflow = parse_project_config(
                {
                    "schema_version": 1,
                    "project": {"id": "empty", "modality": "bulk", "output_root": "results"},
                    "execution": {"profile": "docker", "resume": False},
                    "workflow": {"stages": [{
                        "id": "enrich",
                        "uses": "bulk.enrichment",
                        "inputs": {"de_results": "de.tsv", "genesets": "sets.gmt"},
                    }]},
                },
                source_path=root / "enrichment.toml",
            )
            enrichment = plan_workflow(enrichment_workflow, DEFAULT_CONFIG).stage("enrich")
            enrichment.outdir.mkdir(parents=True)
            self.assertFalse(verify_required_outputs(enrichment.required_outputs).ok)
            (enrichment.outdir / "enrichment-results.tsv").write_text("ok\n", encoding="utf-8")
            self.assertTrue(verify_required_outputs(enrichment.required_outputs).ok)

    def test_rejects_unknown_params_inputs_outputs_and_modality(self) -> None:
        cases = [
            (
                resolved(
                    "bulk",
                    [{"id": "raw", "uses": "bulk.rnaseq", "inputs": {"samplesheet": "s.csv"}, "params": {"typo": True}}],
                ),
                "Unknown parameter",
            ),
            (
                resolved(
                    "bulk",
                    [{"id": "raw", "uses": "bulk.rnaseq", "inputs": {"wrong": "s.csv"}}],
                ),
                "missing required input",
            ),
            (
                resolved(
                    "bulk",
                    [{"id": "raw", "uses": "bulk.rnaseq", "inputs": {"samplesheet": "s.csv"}, "outputs": {"wrong": "x"}}],
                ),
                "unknown output",
            ),
            (
                resolved(
                    "bulk",
                    [{"id": "raw", "uses": "single-cell.scrnaseq", "inputs": {"samplesheet": "s.csv"}}],
                ),
                "modality",
            ),
        ]
        for workflow, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(PlanningError, message):
                plan_workflow(workflow, DEFAULT_CONFIG)

    def test_rejects_unknown_artifacts_and_type_mismatches(self) -> None:
        workflow = ResolvedWorkflow(
            schema_version=1,
            project=ProjectDefinition("x", "single-cell", Path("/work/out")),
            execution=ExecutionSettings("docker", False),
            source_path=SOURCE,
            source_dir=SOURCE.parent,
            stages=(
                StageDefinition(
                    "aggregate",
                    "single-cell.pseudobulk",
                    inputs={"h5ad": Path("/work/a.h5ad")},
                    params={"sample_column": "sample", "design": "group"},
                ),
                StageDefinition(
                    "scan",
                    "single-cell.scanpy",
                    inputs={"input": ArtifactReference("aggregate", "counts")},
                ),
            ),
        )
        with self.assertRaisesRegex(PlanningError, "Artifact type mismatch"):
            plan_workflow(workflow, DEFAULT_CONFIG)

        bad = ResolvedWorkflow(
            schema_version=1,
            project=workflow.project,
            execution=workflow.execution,
            source_path=SOURCE,
            source_dir=SOURCE.parent,
            stages=(
                StageDefinition(
                    "aggregate",
                    "single-cell.pseudobulk",
                    inputs={"h5ad": Path("/work/a.h5ad")},
                    params={"sample_column": "sample", "design": "group"},
                ),
                StageDefinition(
                    "scan",
                    "single-cell.scanpy",
                    inputs={"input": ArtifactReference("aggregate", "missing")},
                ),
            ),
        )
        with self.assertRaisesRegex(PlanningError, "unknown artifact"):
            plan_workflow(bad, DEFAULT_CONFIG)

    def test_detects_cycles_in_manually_constructed_resolved_workflow(self) -> None:
        workflow = ResolvedWorkflow(
            schema_version=1,
            project=ProjectDefinition("x", "bulk", Path("/work/out")),
            execution=ExecutionSettings("docker", False),
            source_path=SOURCE,
            source_dir=SOURCE.parent,
            stages=(
                StageDefinition("a", "bulk.rnaseq", depends_on=("b",), inputs={"samplesheet": Path("/work/s.csv")}),
                StageDefinition("b", "bulk.rnaseq", depends_on=("a",), inputs={"samplesheet": Path("/work/s.csv")}),
            ),
        )
        with self.assertRaisesRegex(PlanningError, "cycle"):
            plan_workflow(workflow, DEFAULT_CONFIG)

    def test_plan_is_immutable_serializable_human_readable_and_side_effect_free(self) -> None:
        workflow = resolved(
            "bulk",
            [{"id": "raw", "uses": "bulk.rnaseq", "inputs": {"samplesheet": "s.csv"}}],
        )
        with patch("subprocess.run", side_effect=AssertionError("must not execute")):
            plan = plan_workflow(workflow, DEFAULT_CONFIG)
        self.assertEqual(json.loads(plan.to_json()), plan.to_dict())
        self.assertIn("maturity", plan.to_json())
        self.assertIn("executables", format_plan_human(plan))
        self.assertEqual(render_plan(plan, format="json"), plan.to_json(indent=2))
        with self.assertRaises(FrozenInstanceError):
            plan.resume = False  # type: ignore[misc]
        with self.assertRaises(TypeError):
            plan.stage("raw").params["new"] = True  # type: ignore[index]


if __name__ == "__main__":
    unittest.main()
