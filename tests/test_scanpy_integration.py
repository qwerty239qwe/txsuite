"""Harmony batch integration through the project workflow API.

Harmony has been available from the CLI and in the single-cell image since the
phase-2 work; what was missing was any way to reach it from a workflow.toml.
These tests cover the stage surface and the wiring, not the algorithm, which the
container smoke in CI already exercises against real data.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG, load_config
from txsuite.project.adapters.single_cell import scanpy_command
from txsuite.project.config import parse_project_config
from txsuite.project.planner import PlanningError, plan_workflow
from txsuite.project.registry import get_stage_spec, validate_stage_parameters
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import INTEGRATION_METHODS, analysis_command


SOURCE = Path("/work/workflow.toml")


class ScanpyIntegrationParameterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = get_stage_spec("single-cell.scanpy")

    def test_defaults_leave_integration_off(self) -> None:
        defaults = validate_stage_parameters(self.spec, {})
        self.assertEqual(defaults["integration"], "none")
        self.assertIsNone(defaults["batch_column"])
        self.assertEqual(defaults["barcode_column"], "barcode")

    def test_registry_and_library_accept_the_same_methods(self) -> None:
        # Two places name the method set; drift between them would let a stage
        # accept a value the container then rejects.
        for method in INTEGRATION_METHODS:
            with self.subTest(method=method):
                params = validate_stage_parameters(self.spec, {"integration": method})
                self.assertEqual(params["integration"], method)
        with self.assertRaises(TxSuiteError):
            validate_stage_parameters(self.spec, {"integration": "harmnoy"})

    def test_column_names_are_validated(self) -> None:
        params = validate_stage_parameters(
            self.spec, {"batch_column": "batch", "barcode_column": "cell_id"}
        )
        self.assertEqual(params["batch_column"], "batch")
        self.assertEqual(params["barcode_column"], "cell_id")
        for bad in ({"batch_column": "not a column"}, {"barcode_column": "1st"}):
            with self.subTest(bad=bad), self.assertRaises(TxSuiteError):
                validate_stage_parameters(self.spec, bad)

    def test_metadata_is_an_optional_input(self) -> None:
        self.assertEqual(
            dict(self.spec.optional_inputs), {"metadata": "single-cell.cell-metadata"}
        )
        self.assertEqual(set(self.spec.inputs), {"input"})


class ScanpyIntegrationCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.matrix = self.root / "matrix.h5ad"
        self.matrix.write_bytes(b"h5ad")
        self.metadata = self.root / "cell-metadata.tsv"
        self.metadata.write_text("barcode\tbatch\nAAAC\tA\n", encoding="utf-8")
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _context(self, params, *, with_metadata=True):
        inputs = {"input": self.matrix}
        if with_metadata:
            inputs["metadata"] = self.metadata
        return {
            "global_config": self.config,
            "resolved_inputs": inputs,
            "params": params,
            "outdir": self.outdir,
        }

    def test_harmony_reaches_the_container_command(self) -> None:
        command = scanpy_command(
            self._context({"integration": "harmony", "batch_column": "batch"})
        )
        self.assertEqual(command[command.index("--integration") + 1], "harmony")
        self.assertEqual(command[command.index("--batch-column") + 1], "batch")
        # The metadata file is mounted read-only and referred to by its
        # in-container path, never the host path.
        self.assertEqual(command[command.index("--metadata") + 1], "/input/metadata.tsv")
        self.assertTrue(
            any(
                f"source={self.metadata.resolve()},target=/input/metadata.tsv,readonly" in item
                for item in command
            )
        )

    def test_adapter_matches_the_library_helper(self) -> None:
        self.assertEqual(
            scanpy_command(
                self._context({"integration": "harmony", "batch_column": "batch"})
            ),
            analysis_command(
                image=self.config["images"]["single_cell_python"],
                input_path=self.matrix,
                outdir=self.outdir,
                min_genes=200,
                min_cells=3,
                max_mito_pct=20.0,
                resolution=1.0,
                metadata=self.metadata,
                barcode_column="barcode",
                batch_column="batch",
                integration="harmony",
            ),
        )

    def test_defaults_produce_an_unintegrated_run_without_metadata(self) -> None:
        command = scanpy_command(self._context({}, with_metadata=False))
        self.assertEqual(command[command.index("--integration") + 1], "none")
        self.assertNotIn("--metadata", command)
        self.assertNotIn("--batch-column", command)

    def test_harmony_without_a_batch_column_is_rejected(self) -> None:
        with self.assertRaisesRegex(TxSuiteError, "requires a batch column"):
            scanpy_command(self._context({"integration": "harmony"}))


class ScanpyIntegrationPlanningTests(unittest.TestCase):
    def _plan(self, stage_inputs, params):
        workflow = parse_project_config(
            {
                "schema_version": 1,
                "project": {
                    "id": "harmony",
                    "modality": "single-cell",
                    "output_root": "results",
                },
                "execution": {"profile": "docker", "resume": True},
                "workflow": {
                    "stages": [
                        {
                            "id": "analyze",
                            "uses": "single-cell.scanpy",
                            "inputs": stage_inputs,
                            "params": params,
                        }
                    ]
                },
            },
            source_path=SOURCE,
        )
        return plan_workflow(workflow, DEFAULT_CONFIG)

    def test_workflow_can_enable_harmony_with_a_metadata_file(self) -> None:
        plan = self._plan(
            {"input": "matrix.h5ad", "metadata": "cell-metadata.tsv"},
            {"integration": "harmony", "batch_column": "batch"},
        )
        command = " ".join(plan.stage("analyze").command)
        self.assertIn("--integration harmony", command)
        self.assertIn("--batch-column batch", command)

    def test_metadata_may_be_omitted_when_obs_already_carries_the_column(self) -> None:
        plan = self._plan(
            {"input": "matrix.h5ad"},
            {"integration": "harmony", "batch_column": "batch"},
        )
        analyze = plan.stage("analyze")
        self.assertEqual(set(analyze.inputs), {"input"})
        self.assertIn("--integration harmony", " ".join(analyze.command))

    def test_unknown_inputs_are_still_rejected(self) -> None:
        with self.assertRaises(PlanningError):
            self._plan({"input": "matrix.h5ad", "cells": "cells.tsv"}, {})


if __name__ == "__main__":
    unittest.main()
