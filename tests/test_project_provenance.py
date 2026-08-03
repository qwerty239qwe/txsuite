from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from txsuite.project.provenance import (
    REDACTED,
    RunBundle,
    atomic_write_json,
    generate_run_id,
    fingerprint_path,
    hash_payload,
    read_json,
    redact_arguments,
)
from txsuite.project.results import (
    StageResult,
    StageState,
    decide_resume,
    verify_required_outputs,
)


class ProjectProvenanceTest(unittest.TestCase):
    def test_directory_fingerprint_is_deterministic_and_content_sensitive(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "input"
            nested = root / "nested"
            nested.mkdir(parents=True)
            (root / "a.txt").write_text("a", encoding="utf-8")
            target = nested / "b.txt"
            target.write_text("b", encoding="utf-8")
            first = fingerprint_path(root)
            second = fingerprint_path(root)
            self.assertEqual(first, second)
            self.assertEqual(first["kind"], "directory")
            self.assertEqual(
                [(item["path"], item["kind"]) for item in first["entries"]],
                [
                    ("a.txt", "file"),
                    ("nested", "directory"),
                    ("nested/b.txt", "file"),
                ],
            )
            target.write_text("changed", encoding="utf-8")
            self.assertNotEqual(first["sha256"], fingerprint_path(root)["sha256"])

    def test_run_id_and_atomic_json(self) -> None:
        first = generate_run_id("analysis")
        second = generate_run_id("analysis")
        self.assertNotEqual(first, second)
        self.assertTrue(first.startswith("analysis-"))
        with self.assertRaises(ValueError):
            generate_run_id("../unsafe")

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nested" / "payload.json"
            atomic_write_json(path, {"version": 1}, overwrite=False)
            self.assertEqual(read_json(path), {"version": 1})
            with self.assertRaises(FileExistsError):
                atomic_write_json(path, {"version": 2}, overwrite=False)
            self.assertEqual(read_json(path), {"version": 1})
            self.assertEqual(list(path.parent.glob("*.tmp")), [])

    def test_bundle_layout_is_pending_and_redacts_persisted_commands(self) -> None:
        command_plan = {
            "stages": [
                {
                    "id": "align",
                    "command": [
                        "tool",
                        "--token",
                        "top-secret",
                        "--password=hunter2",
                        "https://alice:password@example.test/data",
                    ],
                    "env": {"API_TOKEN": "also-secret", "THREADS": "2"},
                }
            ]
        }
        with tempfile.TemporaryDirectory() as directory:
            bundle = RunBundle.create(
                Path(directory),
                run_id="run-fixed",
                workflow={
                    "name": "demo",
                    "params": {"nested": {"api_token": "workflow-secret"}},
                },
                resolved_config={
                    "profile": "local",
                    "credentials": {"password": "config-secret"},
                },
                command_plan=command_plan,
            )
            self.assertEqual(
                {path.name for path in bundle.run_dir.iterdir()},
                {
                    "workflow.json",
                    "resolved-config.json",
                    "run-manifest.json",
                    "command-plan.json",
                    "logs",
                    "stage-results",
                },
            )
            self.assertEqual(list((bundle.run_dir / "logs").iterdir()), [])
            self.assertEqual(list((bundle.run_dir / "stage-results").iterdir()), [])
            manifest = bundle.manifest()
            self.assertEqual(manifest["status"], "pending")
            self.assertEqual(manifest["stage_states"], {"align": "pending"})
            self.assertEqual(
                manifest["hashes"]["command_plan"], hash_payload(command_plan)
            )

            persisted = read_json(bundle.run_dir / "command-plan.json")
            command = persisted["stages"][0]["command"]
            self.assertEqual(command[2], REDACTED)
            self.assertEqual(command[3], f"--password={REDACTED}")
            self.assertNotIn("alice:password", command[4])
            self.assertEqual(persisted["stages"][0]["env"]["API_TOKEN"], REDACTED)
            self.assertEqual(bundle.command_plan(), command_plan)
            workflow = read_json(bundle.run_dir / "workflow.json")
            config = read_json(bundle.run_dir / "resolved-config.json")
            self.assertEqual(
                workflow["params"]["nested"]["api_token"], REDACTED
            )
            self.assertEqual(config["credentials"], REDACTED)
            self.assertNotIn(
                "workflow-secret", (bundle.run_dir / "workflow.json").read_text()
            )
            self.assertNotIn(
                "config-secret", (bundle.run_dir / "resolved-config.json").read_text()
            )
            self.assertEqual(
                manifest["hashes"]["workflow"],
                hash_payload(
                    {
                        "name": "demo",
                        "params": {
                            "nested": {"api_token": "workflow-secret"}
                        },
                    }
                ),
            )

            with self.assertRaisesRegex(ValueError, "Unknown stage state"):
                bundle.set_stage_state("align", "successful")
            self.assertEqual(bundle.manifest()["stage_states"]["align"], "pending")

    def test_invalid_plan_does_not_allocate_a_run_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, "must have an 'id'"):
                RunBundle.create(
                    root,
                    run_id="invalid-run",
                    workflow={},
                    resolved_config={},
                    command_plan={"stages": [{"command": ["true"]}]},
                )
            self.assertFalse((root / "invalid-run").exists())
            with self.assertRaisesRegex(ValueError, "non-empty command list"):
                RunBundle.create(
                    root,
                    run_id="invalid-command-run",
                    workflow={},
                    resolved_config={},
                    command_plan={"stages": [{"id": "stage", "command": "true"}]},
                )
            self.assertFalse((root / "invalid-command-run").exists())

    def test_redaction_handles_split_inline_and_environment_arguments(self) -> None:
        self.assertEqual(
            redact_arguments(
                ["program", "--api-key", "abc", "PASSWORD=def", "--threads=4"]
            ),
            ["program", "--api-key", REDACTED, f"PASSWORD={REDACTED}", "--threads=4"],
        )

    def test_output_verification_and_resume_decision(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "result.txt"
            output.write_text("ok", encoding="utf-8")
            verification = verify_required_outputs(
                [
                    {"path": "result.txt", "kind": "file", "non_empty": True},
                    "missing.txt",
                ],
                base_dir=root,
            )
            self.assertFalse(verification.ok)
            self.assertEqual(verification.missing, (str(root / "missing.txt"),))

            previous = StageResult(
                stage_id="stage",
                state=StageState.COMPLETED,
                input_hash="inputs",
                config_hash="config",
                execution_hash="execution",
                required_outputs=("result.txt",),
            )
            decision = decide_resume(
                previous,
                input_hash="inputs",
                config_hash="config",
                execution_hash="execution",
                required_outputs=["result.txt"],
                base_dir=root,
            )
            self.assertTrue(decision.should_skip)
            changed = decide_resume(
                previous,
                input_hash="inputs-v2",
                config_hash="config",
                execution_hash="execution-v2",
                required_outputs=["result.txt"],
                base_dir=root,
            )
            self.assertTrue(changed.should_run)
            self.assertIn("inputs changed", changed.reason)


if __name__ == "__main__":
    unittest.main()
