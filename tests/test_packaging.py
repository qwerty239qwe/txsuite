from __future__ import annotations

import os
import tomllib
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).parents[1]
PROJECT_RESOURCES = ROOT / "src" / "txsuite" / "resources" / "project"


class PackagingTests(unittest.TestCase):
    def test_wheel_configuration_explicitly_includes_project_resources(self) -> None:
        configuration = tomllib.loads(
            (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        )
        wheel = configuration["tool"]["hatch"]["build"]["targets"]["wheel"]
        self.assertEqual(wheel["packages"], ["src/txsuite"])
        self.assertTrue((PROJECT_RESOURCES / "workflow.schema.json").is_file())
        self.assertTrue((PROJECT_RESOURCES / "presets").is_dir())

    def test_built_wheel_contains_schema_and_every_preset_file(self) -> None:
        wheel_value = os.environ.get("TXSUITE_WHEEL")
        if not wheel_value:
            self.skipTest("set TXSUITE_WHEEL to run the offline wheel audit")
        wheel_path = Path(wheel_value)
        self.assertTrue(wheel_path.is_file(), wheel_path)
        expected = {
            "txsuite/resources/project/"
            + path.relative_to(PROJECT_RESOURCES).as_posix()
            for path in PROJECT_RESOURCES.rglob("*")
            if path.is_file() and "__pycache__" not in path.parts
        }
        with zipfile.ZipFile(wheel_path) as archive:
            members = set(archive.namelist())
        expected.add("txsuite/resources/single_cell_python/collect_de.py")
        self.assertLessEqual(expected, members)


if __name__ == "__main__":
    unittest.main()
