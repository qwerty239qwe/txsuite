from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from txsuite.cli import _effective_project_stages, run
from txsuite.project import read_json, scaffold_project_preset


class _FakeWorkflow:
    def __init__(self, root: Path, project_id: str = "cli_test") -> None:
        self.source_path = root / "workflow.toml"
        self.project = SimpleNamespace(
            id=project_id,
            output_root=root / "outputs",
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "source_path": str(self.source_path),
            "project": {
                "id": self.project.id,
                "output_root": str(self.project.output_root),
            },
        }


class _FakePlan:
    def __init__(self, stages: list[dict[str, object]]) -> None:
        self.stages = tuple(stages)

    def to_dict(self) -> dict[str, object]:
        return {"schema_version": 1, "stages": list(self.stages)}


class ProjectCliTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def invoke(self, *arguments: str) -> tuple[int, str]:
        output = StringIO()
        environment = {"XDG_CONFIG_HOME": str(self.root / "config-home")}
        with patch.dict(os.environ, environment), redirect_stdout(output):
            status = run(list(arguments))
        return status, output.getvalue()

    def scaffold(self, name: str = "scrnaseq-pseudobulk") -> Path:
        return scaffold_project_preset(name, self.root / name)

    def test_init_scaffolds_a_preset_without_overwriting(self) -> None:
        target = self.root / "nested" / "project"

        status, output = self.invoke(
            "project", "init", "--preset", "scrnaseq", str(target)
        )

        self.assertEqual(status, 0)
        self.assertEqual(output.strip(), str(target))
        self.assertTrue((target / "workflow.toml").is_file())
        status, output = self.invoke(
            "project", "init", "--preset", "scrnaseq", str(target)
        )
        self.assertEqual(status, 2)
        self.assertIn("already exists", output)

    def test_validate_crosses_the_planning_boundary(self) -> None:
        project = self.scaffold()
        workflow_path = project / "workflow.toml"

        status, output = self.invoke(
            "project", "validate", str(workflow_path), "--config", str(project / "global.toml")
        )
        self.assertEqual(status, 0)
        self.assertIn("Valid project workflow", output)
        self.assertIn("1 stage(s)", output)

        document = workflow_path.read_text(encoding="utf-8")
        workflow_path.write_text(
            document.replace(
                'uses = "single-cell.pseudobulk-de"',
                'uses = "single-cell.not-registered"',
            ),
            encoding="utf-8",
        )
        status, output = self.invoke("project", "validate", str(workflow_path))
        self.assertEqual(status, 2)
        self.assertIn("Unknown project stage", output)

    def test_plan_json_uses_the_requested_global_config(self) -> None:
        project = self.scaffold()
        config = project / "custom.toml"
        config.write_text(
            '[images]\n'
            'single_cell_python = "registry.example/single-cell:locked"\n'
            'bulk_r = "registry.example/bulk:locked"\n',
            encoding="utf-8",
        )

        status, output = self.invoke(
            "project",
            "plan",
            str(project / "workflow.toml"),
            "--config",
            str(config),
            "--json",
        )

        self.assertEqual(status, 0)
        payload = json.loads(output)
        self.assertEqual(payload["project"]["id"], "scrnaseq_pseudobulk_example")
        self.assertIn(
            "registry.example/single-cell:locked",
            payload["requirements"]["images"],
        )
        self.assertIn("registry.example/bulk:locked", payload["requirements"]["images"])

    def test_invalid_global_config_is_an_actionable_error(self) -> None:
        project = self.scaffold()
        config = project / "invalid.toml"
        config.write_text('[execution]\nprofile = "local"\n', encoding="utf-8")

        status, output = self.invoke(
            "project", "plan", str(project / "workflow.toml"), "--config", str(config)
        )

        self.assertEqual(status, 2)
        self.assertIn("execution.profile", output)

    def test_dry_run_and_invalid_selection_do_not_write(self) -> None:
        project = self.scaffold()
        output_root = project / "results"
        command = (
            "project",
            "run",
            str(project / "workflow.toml"),
            "--config",
            str(project / "global.toml"),
            "--dry-run",
            "--run-id",
            "never-created",
        )

        status, output = self.invoke(*command)

        self.assertEqual(status, 0)
        self.assertIn("Project scrnaseq_pseudobulk_example", output)
        self.assertFalse(output_root.exists())

        status, output = self.invoke(*command, "--stages", "missing")
        self.assertEqual(status, 2)
        self.assertIn("Unknown selected stage", output)
        self.assertFalse(output_root.exists())

    def test_resume_requires_an_existing_named_run_without_writing(self) -> None:
        project = self.scaffold()

        status, output = self.invoke(
            "project", "run", str(project / "workflow.toml"), "--resume"
        )

        self.assertEqual(status, 2)
        self.assertIn("--resume requires --run-id", output)
        self.assertFalse((project / "results").exists())

        status, output = self.invoke(
            "project",
            "run",
            str(project / "workflow.toml"),
            "--resume",
            "--run-id",
            "../outside",
        )
        self.assertEqual(status, 2)
        self.assertIn("Unsafe run ID", output)
        self.assertFalse((project / "results").exists())

    def test_run_status_and_resume_use_the_canonical_bundle(self) -> None:
        workflow = _FakeWorkflow(self.root)
        workspace = self.root / "workspace"
        workspace.mkdir()
        counter = workspace / "counter.txt"
        script = (
            "from pathlib import Path; "
            "p=Path('counter.txt'); "
            "p.write_text(str(int(p.read_text()) + 1) if p.exists() else '1')"
        )
        plan = _FakePlan(
            [
                {
                    "id": "prepare",
                    "command": [sys.executable, "-c", script],
                    "cwd": str(workspace),
                    "required_outputs": [],
                }
            ]
        )
        context = (workflow, {"execution": {"profile": "docker"}}, plan)

        with patch("txsuite.cli._load_project_plan", return_value=context):
            status, output = self.invoke(
                "project", "run", str(workflow.source_path), "--run-id", "fixed"
            )

        run_dir = workflow.project.output_root / workflow.project.id / "runs" / "fixed"
        self.assertEqual(status, 0)
        self.assertEqual(output.splitlines()[0], str(run_dir))
        self.assertIn("Selected stages: prepare", output)
        self.assertEqual(counter.read_text(encoding="utf-8"), "1")
        self.assertEqual(read_json(run_dir / "run-manifest.json")["status"], "completed")

        manifest_before = (run_dir / "run-manifest.json").read_bytes()
        status, output = self.invoke("project", "status", str(run_dir))
        self.assertEqual(status, 0)
        self.assertEqual(json.loads(output)["run_id"], "fixed")
        self.assertEqual((run_dir / "run-manifest.json").read_bytes(), manifest_before)

        with patch("txsuite.cli._load_project_plan", return_value=context):
            status, output = self.invoke(
                "project",
                "run",
                str(workflow.source_path),
                "--resume",
                "--run-id",
                "fixed",
            )
        self.assertEqual(status, 0)
        self.assertEqual(counter.read_text(encoding="utf-8"), "1")
        result = read_json(run_dir / "stage-results" / "prepare.json")
        self.assertEqual(len(result["attempts"]), 1)
        self.assertTrue(result["reason"].startswith("resume:"))

    def test_stage_selection_and_failure_are_recorded(self) -> None:
        workflow = _FakeWorkflow(self.root, "selection")
        workspace = self.root / "selection-workspace"
        workspace.mkdir()
        plan = _FakePlan(
            [
                {
                    "id": "first",
                    "command": [sys.executable, "-c", "open('first', 'w').close()"],
                    "cwd": str(workspace),
                    "required_outputs": [],
                },
                {
                    "id": "bad",
                    "command": [sys.executable, "-c", "import sys; sys.exit(7)"],
                    "cwd": str(workspace),
                    "required_outputs": [],
                },
            ]
        )
        context = (workflow, {}, plan)

        with patch("txsuite.cli._load_project_plan", return_value=context):
            status, output = self.invoke(
                "project",
                "run",
                str(workflow.source_path),
                "--run-id",
                "failed-run",
                "--stages",
                "bad",
            )

        run_dir = workflow.project.output_root / "selection" / "runs" / "failed-run"
        self.assertEqual(status, 2)
        self.assertIn("exited with status 7", output)
        self.assertIn(str(run_dir), output)
        self.assertFalse((workspace / "first").exists())
        manifest = read_json(run_dir / "run-manifest.json")
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(manifest["stage_states"]["first"], "skipped")
        self.assertEqual(manifest["stage_states"]["bad"], "failed")

    def test_dependency_closure_preserves_order_and_excludes_other_branches(self) -> None:
        plan = _FakePlan(
            [
                {"id": "raw", "command": ["true"], "depends_on": []},
                {"id": "independent", "command": ["true"], "depends_on": []},
                {"id": "de", "command": ["true"], "depends_on": ["raw"]},
                {"id": "report", "command": ["true"], "depends_on": ["de"]},
            ]
        )

        effective, auto_added = _effective_project_stages(
            plan,
            selected_stage_ids=("de",),
            from_stage=None,
            to_stage=None,
        )
        self.assertEqual(effective, ("raw", "de"))
        self.assertEqual(auto_added, ("raw",))

        effective, auto_added = _effective_project_stages(
            plan,
            selected_stage_ids=None,
            from_stage="de",
            to_stage="report",
        )
        self.assertEqual(effective, ("raw", "de", "report"))
        self.assertEqual(auto_added, ("raw",))
        self.assertNotIn("independent", effective)

    def test_dry_run_shows_dependency_closure_and_omits_unselected_stage(self) -> None:
        project = self.scaffold("bulk-rnaseq")
        output_root = project / "results"

        status, output = self.invoke(
            "project",
            "run",
            str(project / "workflow.toml"),
            "--config",
            str(project / "global.toml"),
            "--dry-run",
            "--stages",
            "differential",
        )

        self.assertEqual(status, 0)
        self.assertIn("Selected stages: rnaseq, differential", output)
        self.assertIn("Auto-added prerequisites: rnaseq", output)
        self.assertIn("bulk.rnaseq", output)
        self.assertIn("bulk.de", output)
        self.assertNotIn("bulk.enrichment", output)
        self.assertFalse(output_root.exists())

    def test_executor_receives_expanded_dependency_ids(self) -> None:
        workflow = _FakeWorkflow(self.root, "expanded")
        plan = _FakePlan(
            [
                {"id": "raw", "command": ["true"], "depends_on": []},
                {"id": "branch", "command": ["true"], "depends_on": []},
                {"id": "de", "command": ["true"], "depends_on": ["raw"]},
            ]
        )
        context = (workflow, {}, plan)

        with (
            patch("txsuite.cli._load_project_plan", return_value=context),
            patch("txsuite.cli.ProjectExecutor.execute") as execute,
        ):
            status, output = self.invoke(
                "project",
                "run",
                str(workflow.source_path),
                "--run-id",
                "selected",
                "--stages",
                "de",
            )

        self.assertEqual(status, 0)
        self.assertIn("Selected stages: raw, de", output)
        self.assertIn("Auto-added prerequisites: raw", output)
        execute.assert_called_once()
        self.assertEqual(
            execute.call_args.kwargs["selected_stage_ids"], ("raw", "de")
        )
        self.assertNotIn("from_stage", execute.call_args.kwargs)
        self.assertNotIn("to_stage", execute.call_args.kwargs)

    def test_resume_rejects_changed_snapshots_without_execution_or_mutation(self) -> None:
        workflow = _FakeWorkflow(self.root, "immutable")
        workspace = self.root / "immutable-workspace"
        workspace.mkdir()
        original_plan = _FakePlan(
            [
                {
                    "id": "stage",
                    "command": [sys.executable, "-c", "print('original')"],
                    "cwd": str(workspace),
                    "required_outputs": [],
                }
            ]
        )
        original_context = (
            workflow,
            {"execution": {"profile": "docker"}},
            original_plan,
        )
        with patch("txsuite.cli._load_project_plan", return_value=original_context):
            status, _ = self.invoke(
                "project", "run", str(workflow.source_path), "--run-id", "fixed"
            )
        self.assertEqual(status, 0)

        run_dir = workflow.project.output_root / "immutable" / "runs" / "fixed"
        before = {
            path.relative_to(run_dir).as_posix(): path.read_bytes()
            for path in run_dir.rglob("*")
            if path.is_file()
        }
        changed_plan = _FakePlan(
            [
                {
                    "id": "stage",
                    "command": [
                        sys.executable,
                        "-c",
                        "open('must-not-run', 'w').close()",
                    ],
                    "cwd": str(workspace),
                    "required_outputs": [],
                }
            ]
        )
        changed_context = (
            workflow,
            {"execution": {"profile": "apptainer"}},
            changed_plan,
        )
        with (
            patch("txsuite.cli._load_project_plan", return_value=changed_context),
            patch("txsuite.project.executor.subprocess.run") as subprocess_run,
        ):
            status, output = self.invoke(
                "project",
                "run",
                str(workflow.source_path),
                "--resume",
                "--run-id",
                "fixed",
            )

        self.assertEqual(status, 2)
        self.assertIn("resolved configuration", output)
        self.assertIn("command plan", output)
        self.assertIn("Start a new run", output)
        subprocess_run.assert_not_called()
        self.assertFalse((workspace / "must-not-run").exists())
        after = {
            path.relative_to(run_dir).as_posix(): path.read_bytes()
            for path in run_dir.rglob("*")
            if path.is_file()
        }
        self.assertEqual(after, before)

    def test_status_normalizes_missing_and_invalid_bundles(self) -> None:
        missing = self.root / "missing-run"
        status, output = self.invoke("project", "status", str(missing))
        self.assertEqual(status, 2)
        self.assertIn("Cannot read run bundle", output)

        invalid = self.root / "invalid-run"
        invalid.mkdir()
        (invalid / "run-manifest.json").write_text("not JSON\n", encoding="utf-8")
        status, output = self.invoke("project", "status", str(invalid))
        self.assertEqual(status, 2)
        self.assertIn("Cannot read run bundle", output)


if __name__ == "__main__":
    unittest.main()
