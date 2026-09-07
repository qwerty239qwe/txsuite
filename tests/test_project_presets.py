from __future__ import annotations

import os
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG
from txsuite.project.config import load_project_config
from txsuite.project.planner import plan_workflow
from txsuite.project.presets import (
    get_project_preset,
    list_project_presets,
    scaffold_project_preset,
)
from txsuite.runtime import TxSuiteError


EXPECTED_PRESETS = (
    "bulk-rnaseq",
    "scrnaseq",
    "scrnaseq-pseudobulk",
)


class ProjectPresetTests(unittest.TestCase):
    def test_lists_and_gets_immutable_presets_in_stable_order(self) -> None:
        presets = list_project_presets()

        self.assertEqual(tuple(preset.name for preset in presets), EXPECTED_PRESETS)
        self.assertEqual(get_project_preset("bulk-rnaseq"), presets[0])
        self.assertEqual(presets[0].id, "bulk-rnaseq")
        self.assertEqual(presets[0].modality, "bulk")
        with self.assertRaises(FrozenInstanceError):
            presets[0].name = "changed"  # type: ignore[misc]

    def test_unknown_preset_reports_the_available_names(self) -> None:
        with self.assertRaisesRegex(TxSuiteError, "Unknown project preset") as raised:
            get_project_preset("unknown")
        for name in EXPECTED_PRESETS:
            self.assertIn(name, str(raised.exception))

    def test_scaffolds_each_preset_and_can_load_and_plan_it(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            unrelated_cwd = root / "unrelated-cwd"
            unrelated_cwd.mkdir()
            previous_cwd = Path.cwd()
            try:
                os.chdir(unrelated_cwd)
                for preset in list_project_presets():
                    with self.subTest(preset=preset.name):
                        target = root / "nested" / preset.name
                        result = scaffold_project_preset(preset.name, target)

                        self.assertEqual(result, target)
                        self.assertEqual(
                            tuple(
                                path.relative_to(target).as_posix()
                                for path in sorted(target.rglob("*"), key=lambda p: p.relative_to(target).as_posix())
                                if path.is_file()
                            ),
                            preset.files,
                        )
                        workflow = load_project_config(target / "workflow.toml")
                        self.assertEqual(workflow.project.modality, preset.modality)

                        plan = plan_workflow(workflow, DEFAULT_CONFIG)
                        self.assertEqual(plan.project_id, workflow.project.id)
                        self.assertEqual(
                            tuple(stage.id for stage in plan.stages),
                            tuple(stage.id for stage in workflow.stages),
                        )
                        params_files = {
                            "bulk-rnaseq": "rnaseq-params.json",
                            "scrnaseq": "scrnaseq-params.json",
                        }
                        if preset.name in params_files:
                            command = plan.stages[0].command
                            index = command.index("-params-file")
                            self.assertEqual(
                                command[index + 1],
                                str((target / params_files[preset.name]).resolve()),
                            )
            finally:
                os.chdir(previous_cwd)

    def test_existing_target_is_rejected_without_changing_it(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            occupied = root / "occupied"
            occupied.mkdir()
            marker = occupied / "keep.txt"
            marker.write_text("user content\n", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "already exists"):
                scaffold_project_preset("bulk-rnaseq", occupied)

            self.assertEqual(marker.read_text(encoding="utf-8"), "user content\n")
            self.assertEqual(tuple(occupied.iterdir()), (marker,))

    def test_existing_empty_target_is_also_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "empty"
            target.mkdir()

            with self.assertRaises(FileExistsError):
                scaffold_project_preset("scrnaseq", target)
            self.assertEqual(tuple(target.iterdir()), ())


if __name__ == "__main__":
    unittest.main()
