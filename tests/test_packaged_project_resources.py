from __future__ import annotations

import csv
import json
import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.project.config import load_project_config
from txsuite.project.model import ArtifactReference
from txsuite.project.presets import list_project_presets, scaffold_project_preset
from txsuite.project.registry import get_stage_spec, validate_stage_parameters


class PackagedProjectResourceTests(unittest.TestCase):
    def test_resource_manifests_exist_and_match_user_facing_examples(self) -> None:
        package_root = resources.files("txsuite.resources.project").joinpath("presets")
        examples_root = Path(__file__).parents[1] / "examples" / "projects"

        for preset in list_project_presets():
            with self.subTest(preset=preset.name):
                resource_root = package_root.joinpath(preset.name)
                self.assertTrue(resource_root.is_dir())
                self.assertEqual(
                    tuple(
                        child.name
                        for child in sorted(
                            resource_root.iterdir(), key=lambda item: item.name
                        )
                        if child.is_file()
                    ),
                    preset.files,
                )
                for relative in preset.files:
                    packaged = resource_root.joinpath(relative).read_bytes()
                    example = (examples_root / preset.name / relative).read_bytes()
                    self.assertEqual(packaged, example, relative)

    def test_workflows_use_registry_contracts_and_valid_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for preset in list_project_presets():
                with self.subTest(preset=preset.name):
                    target = scaffold_project_preset(preset.name, root / preset.name)
                    workflow = load_project_config(target / "workflow.toml")
                    self.assertEqual(workflow.project.modality, preset.modality)

                    for stage in workflow.stages:
                        spec = get_stage_spec(stage.uses)
                        self.assertEqual(spec.modality, workflow.project.modality)
                        # Every required input must be present, and nothing
                        # beyond the optional ones the stage accepts.
                        self.assertLessEqual(set(spec.inputs), set(stage.inputs))
                        self.assertLessEqual(
                            set(stage.inputs), set(spec.accepted_inputs)
                        )
                        self.assertLessEqual(set(stage.outputs), set(spec.outputs))
                        validate_stage_parameters(spec, stage.params)

                        for input_name, value in stage.inputs.items():
                            if not isinstance(value, ArtifactReference):
                                continue
                            producer = workflow.stage(value.stage_id)
                            producer_spec = get_stage_spec(producer.uses)
                            self.assertIn(value.artifact_id, producer_spec.outputs)
                            self.assertEqual(
                                producer_spec.outputs[value.artifact_id],
                                spec.accepted_inputs[input_name],
                            )

    def test_readmes_contain_exact_validate_plan_and_run_commands(self) -> None:
        package_root = resources.files("txsuite.resources.project").joinpath("presets")
        commands = (
            "txsuite project validate workflow.toml",
            "txsuite project plan workflow.toml",
            "txsuite project run workflow.toml",
        )
        for preset in list_project_presets():
            readme = package_root.joinpath(preset.name, "README.md").read_text(
                encoding="utf-8"
            )
            with self.subTest(preset=preset.name):
                self.assertIn("Replace", readme)
                for command in commands:
                    self.assertEqual(readme.count(command), 1)

    def test_tabular_and_json_templates_are_syntactically_valid(self) -> None:
        package_root = resources.files("txsuite.resources.project").joinpath("presets")
        for preset in list_project_presets():
            root = package_root.joinpath(preset.name)
            with self.subTest(preset=preset.name):
                if "samplesheet.csv" in preset.files:
                    samples = list(
                        csv.reader(root.joinpath("samplesheet.csv").read_text(
                            encoding="utf-8"
                        ).splitlines())
                    )
                    self.assertGreaterEqual(len(samples), 2)
                    self.assertTrue(all(cell for cell in samples[0]))
                    self.assertTrue(
                        all(len(row) == len(samples[0]) for row in samples)
                    )

                for relative in preset.files:
                    if relative.endswith("-params.json"):
                        value = json.loads(root.joinpath(relative).read_text("utf-8"))
                        self.assertIsInstance(value, dict)


if __name__ == "__main__":
    unittest.main()
