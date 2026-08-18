"""Optional stage inputs.

Schema v1 originally required every input a stage declared, which made
"reuse a prebuilt artifact when one is supplied" unrepresentable: as a
required input it forces every project to produce one, and as a parameter it
loses ``${stage.artifact}`` resolution, since references are only parsed
inside ``inputs``. Optional inputs close that gap without weakening the
strictness that applies to required inputs and unknown keys.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from txsuite.project.config import parse_project_config
from txsuite.project.planner import PlanningError, plan_workflow
from txsuite.project.registry import StageSpec, get_stage_spec

from txsuite.config import DEFAULT_CONFIG


SOURCE = Path("/work/workflow.toml")


def _spec(**overrides) -> StageSpec:
    base = {
        "id": "bulk.example",
        "modality": "bulk",
        "maturity": "ready",
        "inputs": {"counts": "bulk.gene-counts"},
        "outputs": {"results": "bulk.differential-expression-results"},
        "defaults": {},
        "validators": {},
        "command_factory": lambda context: ["true"],
    }
    base.update(overrides)
    return StageSpec(**base)


class OptionalInputContractTests(unittest.TestCase):
    def test_accepted_inputs_merges_required_and_optional(self) -> None:
        spec = _spec(optional_inputs={"index": "bulk.star-index"})
        self.assertEqual(set(spec.inputs), {"counts"})
        self.assertEqual(set(spec.optional_inputs), {"index"})
        self.assertEqual(set(spec.accepted_inputs), {"counts", "index"})

    def test_an_input_cannot_be_required_and_optional(self) -> None:
        with self.assertRaises(ValueError):
            _spec(optional_inputs={"counts": "bulk.gene-counts"})

    def test_optional_inputs_default_to_empty_and_stay_frozen(self) -> None:
        spec = _spec()
        self.assertEqual(dict(spec.optional_inputs), {})
        with self.assertRaises(TypeError):
            spec.optional_inputs["index"] = "bulk.star-index"  # type: ignore[index]

    def test_optional_contracts_reject_globs_like_required_ones(self) -> None:
        with self.assertRaises(ValueError):
            _spec(optional_inputs={"index": "bulk.*-index"})

    def test_existing_stages_declare_no_optional_inputs(self) -> None:
        for uses in ("bulk.rnaseq", "bulk.de", "single-cell.scanpy"):
            with self.subTest(uses=uses):
                spec = get_stage_spec(uses)
                self.assertEqual(dict(spec.optional_inputs), {})
                self.assertEqual(
                    dict(spec.accepted_inputs), dict(spec.inputs)
                )


class OptionalInputPlanningTests(unittest.TestCase):
    """Plan real workflows through a stage that declares an optional input."""

    def _plan(self, stages):
        workflow = parse_project_config(
            {
                "schema_version": 1,
                "project": {
                    "id": "optional_inputs",
                    "modality": "bulk",
                    "output_root": "results",
                },
                "execution": {"profile": "docker", "resume": True},
                "workflow": {"stages": stages},
            },
            source_path=SOURCE,
        )
        return plan_workflow(workflow, DEFAULT_CONFIG)

    def _stages(self, variant_inputs=None, *, with_reference=False):
        """A bulk.rnavar stage, the first stage to declare optional inputs."""

        inputs = {"samplesheet": "samplesheet.csv"}
        inputs.update(variant_inputs or {})
        stages = [
            {
                "id": "variants",
                "uses": "bulk.rnavar",
                "inputs": inputs,
                "params": {
                    "fasta": "genome.fa",
                    "gtf": "genes.gtf",
                    "skip_baserecalibration": True,
                },
            }
        ]
        if with_reference:
            stages.insert(
                0,
                {
                    "id": "ref",
                    "uses": "bulk.star-reference",
                    "inputs": {"fasta": "genome.fa", "gtf": "genes.gtf"},
                },
            )
        return stages

    def test_planning_succeeds_when_optional_inputs_are_omitted(self) -> None:
        plan = self._plan(self._stages())
        self.assertEqual(set(plan.stage("variants").inputs), {"samplesheet"})
        self.assertNotIn("--star_index", plan.stage("variants").command)

    def test_supplied_optional_inputs_resolve_references_and_imply_an_edge(self) -> None:
        plan = self._plan(
            self._stages(
                {
                    "star_index": "${ref.star_index}",
                    "fasta_fai": "${ref.fasta_fai}",
                    "dict": "${ref.dict}",
                },
                with_reference=True,
            )
        )

        variants = plan.stage("variants")
        reference = plan.stage("ref")
        self.assertIn("ref", variants.depends_on)
        self.assertEqual(
            variants.inputs["star_index"].value, reference.outputs["star_index"].path
        )
        command = " ".join(variants.command)
        self.assertIn(str(reference.outputs["star_index"].path), command)
        self.assertIn(str(reference.outputs["dict"].path), command)

    def test_optional_inputs_are_type_checked_like_required_ones(self) -> None:
        with self.assertRaises(PlanningError) as caught:
            self._plan(
                self._stages({"star_index": "${ref.dict}"}, with_reference=True)
            )
        self.assertIn("Artifact type mismatch", str(caught.exception))

    def test_unknown_inputs_are_still_rejected(self) -> None:
        with self.assertRaises(PlanningError) as caught:
            self._plan(self._stages({"nonsense": "somewhere.tsv"}))
        self.assertIn("unknown input", str(caught.exception))

    def test_missing_required_inputs_are_still_rejected(self) -> None:
        stages = [
            {
                "id": "variants",
                "uses": "bulk.rnavar",
                "inputs": {"star_index": "index"},
                "params": {
                    "fasta": "genome.fa",
                    "gtf": "genes.gtf",
                    "skip_baserecalibration": True,
                },
            }
        ]
        with self.assertRaises(PlanningError) as caught:
            self._plan(stages)
        self.assertIn("missing required input", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
