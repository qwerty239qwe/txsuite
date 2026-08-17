from __future__ import annotations

import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from txsuite.catalog import select_tools
from txsuite.cli import run
from txsuite.config import load_config


class FoundationTest(unittest.TestCase):
    def test_catalog_filter_and_project_config_override(self) -> None:
        self.assertEqual(
            [tool.name for tool in select_tools("bulk", "alignment")], ["STAR"]
        )
        with tempfile.TemporaryDirectory() as directory:
            config_path = Path(directory) / "txsuite.toml"
            config_path.write_text(
                '[execution]\nprofile = "apptainer"\n', encoding="utf-8"
            )
            config = load_config(
                config_path, user_path=Path(directory) / "missing.toml"
            )
        self.assertEqual(config["execution"]["profile"], "apptainer")
        self.assertEqual(config["pipelines"]["bulk"]["release"], "3.26.0")

    def test_env_list_reports_configured_owned_images(self) -> None:
        output = StringIO()
        with redirect_stdout(output):
            status = run(["env", "list", "--config", "txsuite.toml"])

        self.assertEqual(status, 0)
        self.assertEqual(
            output.getvalue().splitlines(),
            [
                "bulk-r\ttxsuite/bulk-r:0.2.0",
                "salmon\ttxsuite/salmon:0.2.0",
                "single-cell-python\ttxsuite/single-cell-python:0.2.0",
                "spatial-python\ttxsuite/spatial-python:0.2.0",
            ],
        )


if __name__ == "__main__":
    unittest.main()
