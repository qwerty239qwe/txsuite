from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG
from txsuite.project.config import parse_project_config
from txsuite.project.executor import (
    ProjectExecutor,
    StageExecutionError,
    select_stages,
)
from txsuite.project.provenance import REDACTED, RunBundle, read_json
from txsuite.project.planner import PlannedWorkflow, plan_workflow
from txsuite.project.results import StageState


def python_stage(
    stage_id: str,
    script: str,
    *,
    cwd: Path | None = None,
    outputs: list[object] | None = None,
    extra_arguments: list[str] | None = None,
    inputs: dict[str, object] | None = None,
) -> dict[str, object]:
    return {
        "id": stage_id,
        "command": [sys.executable, "-c", script, *(extra_arguments or [])],
        "cwd": str(cwd) if cwd is not None else None,
        "required_outputs": outputs or [],
        "inputs": inputs or {},
    }


class ProjectExecutorTest(unittest.TestCase):
    def test_file_parameter_content_affects_resume_fingerprint(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metadata = root / "metadata.tsv"
            metadata.write_text("first")
            executor = ProjectExecutor(self.make_bundle(root, {"stages": []}))
            stage = {"params": {"metadata": str(metadata)}}
            before = executor._input_fingerprints(stage, {}, cwd=root)
            metadata.write_text("other")
            after = executor._input_fingerprints(stage, {}, cwd=root)
            self.assertNotEqual(before, after)

    def make_bundle(
        self, root: Path, command_plan: object, run_id: str = "test-run"
    ) -> RunBundle:
        return RunBundle.create(
            root / "runs",
            run_id=run_id,
            workflow={"name": "test-workflow"},
            resolved_config={"profile": "test"},
            command_plan=command_plan,
        )

    def test_selection_is_deterministic_and_validates_ranges(self) -> None:
        stages = [
            {"id": "one", "command": ["true"]},
            {"id": "two", "command": ["true"]},
            {"id": "three", "command": ["true"]},
        ]
        self.assertEqual(
            [
                stage["id"]
                for stage in select_stages(
                    stages, selected_stage_ids=["three", "one"]
                )
            ],
            ["one", "three"],
        )
        self.assertEqual(
            [
                stage["id"]
                for stage in select_stages(
                    stages, from_stage="two", to_stage="three"
                )
            ],
            ["two", "three"],
        )
        with self.assertRaisesRegex(ValueError, "must not follow"):
            select_stages(stages, from_stage="three", to_stage="one")
        with self.assertRaisesRegex(ValueError, "Unknown selected"):
            select_stages(stages, selected_stage_ids=["missing"])

    def test_accepts_planned_workflow_objects(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = PlannedWorkflow(
                schema_version=1,
                project_id="planned",
                modality="bulk",
                output_root=root / "output",
                profile="local",
                resume=False,
                source_path=root / "workflow.toml",
                stages=(),
            )
            bundle = self.make_bundle(root, plan)
            self.assertEqual(bundle.command_plan(), plan.to_dict())
            self.assertTrue(ProjectExecutor(bundle).execute(plan).succeeded)

    def test_scanpy_symbolic_h5ad_mount_materializes_with_stable_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workflow = parse_project_config(
                {
                    "schema_version": 1,
                    "project": {
                        "id": "scrna",
                        "modality": "single-cell",
                        "output_root": "results",
                    },
                    "execution": {"profile": "docker", "resume": False},
                    "workflow": {
                        "stages": [
                            {
                                "id": "quantify",
                                "uses": "single-cell.scrnaseq",
                                "inputs": {"samplesheet": "samples.csv"},
                            },
                            {
                                "id": "scanpy",
                                "uses": "single-cell.scanpy",
                                "inputs": {"input": "${quantify.matrix}"},
                            },
                        ]
                    },
                },
                source_path=root / "workflow.toml",
            )
            plan = plan_workflow(workflow, DEFAULT_CONFIG)
            preview = plan.stage("scanpy")
            self.assertIn(
                "type=bind,source=${quantify.matrix},target=/input/data.h5ad,readonly",
                preview.command,
            )

            combined = (
                root
                / "results"
                / "quantify"
                / "simpleaf"
                / "mtx_conversions"
                / "combined_matrix.h5ad"
            ).resolve()
            bundle = self.make_bundle(root, plan)
            materialized = ProjectExecutor(bundle)._materialize_stage(
                preview.to_dict(), {"quantify.matrix": {"path": str(combined)}}
            )
            command = materialized["command"]
            self.assertIn(
                f"type=bind,source={combined},target=/input/data.h5ad,readonly",
                command,
            )
            analyze = command.index("analyze")
            self.assertEqual(command[analyze + 1], "/input/data.h5ad")

    def test_success_captures_logs_timing_outputs_and_redacted_command(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workspace = root / "work"
            workspace.mkdir()
            plan = {
                "stages": [
                    python_stage(
                        "prepare",
                        "import pathlib,sys; "
                        "pathlib.Path('artifact.txt').write_text('data'); "
                        "print('hello'); print('warning', file=sys.stderr)",
                        cwd=workspace,
                        outputs=[
                            {"path": "artifact.txt", "kind": "file", "non_empty": True}
                        ],
                        extra_arguments=["--token", "super-secret"],
                        inputs={"sample": "A"},
                    )
                ]
            }
            bundle = self.make_bundle(root, plan)
            result = ProjectExecutor(bundle).execute()

            self.assertTrue(result.succeeded)
            stage = result.stages[0]
            self.assertEqual(stage.state, StageState.COMPLETED)
            self.assertEqual(stage.attempts[0].exit_code, 0)
            self.assertGreaterEqual(stage.attempts[0].duration_seconds or -1, 0)
            self.assertEqual(stage.attempts[0].command[-1], REDACTED)
            self.assertEqual(
                (bundle.run_dir / stage.attempts[0].stdout_log).read_text().strip(),
                "hello",
            )
            self.assertEqual(
                (bundle.run_dir / stage.attempts[0].stderr_log).read_text().strip(),
                "warning",
            )
            self.assertEqual(bundle.manifest()["status"], "completed")
            self.assertEqual(bundle.manifest()["stage_states"]["prepare"], "completed")

            persisted = read_json(bundle.stage_result_path("prepare"))
            serialized = json.dumps(persisted)
            self.assertNotIn("super-secret", serialized)
            self.assertIn(REDACTED, serialized)

    def test_exit_failure_is_recorded_before_raise_and_blocks_following(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = {
                "stages": [
                    python_stage("bad", "import sys; print('bad'); sys.exit(7)"),
                    python_stage("later", "print('should not run')"),
                ]
            }
            bundle = self.make_bundle(root, plan)
            result = ProjectExecutor(bundle).execute()
            self.assertEqual(result.status, "failed")
            self.assertEqual(result.stages[0].state, StageState.FAILED)
            self.assertEqual(result.stages[0].attempts[0].exit_code, 7)
            self.assertEqual(result.stages[1].state, StageState.SKIPPED)
            self.assertIn("blocked by failed stage", result.stages[1].reason or "")
            self.assertEqual(bundle.manifest()["status"], "failed")

            second = self.make_bundle(root, plan, run_id="raising-run")
            with self.assertRaises(StageExecutionError):
                ProjectExecutor(second).execute(raise_on_failure=True)
            self.assertEqual(
                read_json(second.stage_result_path("bad"))["state"], "failed"
            )
            self.assertEqual(second.manifest()["status"], "failed")

    def test_required_output_failure_and_attempt_history(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = {
                "stages": [
                    python_stage(
                        "missing-output",
                        "print('command succeeded')",
                        outputs=["never-created.txt"],
                    )
                ]
            }
            bundle = self.make_bundle(root, plan)
            executor = ProjectExecutor(bundle)
            first = executor.execute()
            self.assertEqual(first.stages[0].state, StageState.FAILED)
            self.assertEqual(first.stages[0].attempts[0].exit_code, 0)
            self.assertIn("required outputs", first.stages[0].reason or "")
            first_record = bundle.attempt_result_path("missing-output", 1)
            self.assertTrue(first_record.is_file())
            first_record_bytes = first_record.read_bytes()
            self.assertIn("missing-output.attempt-1.test-run", first_record.name)
            self.assertEqual(
                first.stages[0].attempt_records,
                (str(first_record.relative_to(bundle.run_dir)),),
            )

            second = executor.execute(resume=True)
            self.assertEqual(second.stages[0].state, StageState.FAILED)
            self.assertEqual(len(second.stages[0].attempts), 2)
            self.assertEqual(first_record.read_bytes(), first_record_bytes)
            second_record = bundle.attempt_result_path("missing-output", 2)
            self.assertTrue(second_record.is_file())
            self.assertNotEqual(first_record, second_record)
            self.assertEqual(len(second.stages[0].attempt_records), 2)
            summary = read_json(bundle.stage_result_path("missing-output"))
            self.assertEqual(
                summary["attempt_records"],
                list(second.stages[0].attempt_records),
            )
            self.assertEqual(
                bundle.manifest()["stage_attempt_records"]["missing-output"],
                list(second.stages[0].attempt_records),
            )

    def test_resume_invalidates_changed_command_or_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workspace = root / "work"
            workspace.mkdir()
            script = (
                "from pathlib import Path; "
                "p=Path('count.txt'); "
                "p.write_text(str(int(p.read_text())+1) if p.exists() else '1')"
            )
            first_stage = python_stage(
                "count",
                script,
                cwd=workspace,
                outputs=["count.txt"],
                inputs={"sample": "A"},
            )
            plan = {"stages": [first_stage]}
            bundle = self.make_bundle(root, plan)
            executor = ProjectExecutor(bundle)
            first = executor.execute(input_hashes={"count": "input-v1"})
            self.assertEqual(first.stages[0].state, StageState.COMPLETED)
            self.assertEqual((workspace / "count.txt").read_text(), "1")

            skipped = executor.execute(resume=True, input_hashes={"count": "input-v1"})
            self.assertEqual(skipped.stages[0].state, StageState.SKIPPED)
            self.assertTrue((skipped.stages[0].reason or "").startswith("resume:"))
            self.assertEqual(len(skipped.stages[0].attempts), 1)
            self.assertEqual((workspace / "count.txt").read_text(), "1")

            changed_command_stage = dict(first_stage)
            changed_command_stage["params"] = {"mode": "changed"}
            changed_plan = {"stages": [changed_command_stage]}
            rerun = executor.execute(
                changed_plan, resume=True, input_hashes={"count": "input-v1"}
            )
            self.assertEqual(rerun.stages[0].state, StageState.COMPLETED)
            self.assertEqual(len(rerun.stages[0].attempts), 2)
            self.assertEqual((workspace / "count.txt").read_text(), "2")

            input_rerun = executor.execute(
                changed_plan, resume=True, input_hashes={"count": "input-v2"}
            )
            self.assertEqual(input_rerun.stages[0].state, StageState.COMPLETED)
            self.assertEqual(len(input_rerun.stages[0].attempts), 3)
            self.assertEqual((workspace / "count.txt").read_text(), "3")

    def test_selected_stage_execution_marks_only_new_unselected_results(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = {
                "stages": [
                    python_stage("one", "print('one')"),
                    python_stage("two", "print('two')"),
                    python_stage("three", "print('three')"),
                ]
            }
            bundle = self.make_bundle(root, plan)
            result = ProjectExecutor(bundle).execute(
                selected_stage_ids=["three", "one"]
            )
            self.assertEqual(
                [stage.state for stage in result.stages],
                [StageState.COMPLETED, StageState.SKIPPED, StageState.COMPLETED],
            )
            self.assertEqual(result.stages[1].reason, "not selected")

    def nfcore_handoff_plan(
        self,
        workspace: Path,
        *,
        mode: str = "unique",
        override: str | None = None,
    ) -> tuple[dict[str, object], Path]:
        results_root = workspace / "rnaseq-results"
        downstream_output = workspace / "downstream.txt"
        raw_script = (
            "from pathlib import Path; import sys; "
            "root=Path(sys.argv[1]); root.mkdir(parents=True, exist_ok=True); "
            "paths=[] if sys.argv[2]=='missing' else "
            "[root/'star_salmon'/'salmon.merged.gene_counts.tsv']; "
            "paths += [root/'salmon'/'salmon.merged.gene_counts.tsv'] "
            "if sys.argv[2]=='ambiguous' else []; "
            "[(p.parent.mkdir(parents=True, exist_ok=True), "
            "p.write_text('' if sys.argv[2]=='zero' else "
            "'gene\\tcount\\nA\\t1\\n')) for p in paths]"
        )
        postflight: dict[str, object] = {
            "adapter": "nfcore",
            "pipeline": "nf-core/rnaseq",
            "release": "3.26.0",
            "results_root": str(results_root),
            "artifacts": {"counts": "bulk.gene-counts"},
            "policies": {"counts": {"kind": "file", "non_empty": True}},
            "overrides": {},
        }
        if override is not None:
            postflight["overrides"] = {"counts": override}
        plan: dict[str, object] = {
            "stages": [
                {
                    "id": "raw",
                    "command": [
                        sys.executable,
                        "-c",
                        raw_script,
                        str(results_root),
                        mode,
                    ],
                    "cwd": str(workspace),
                    "required_outputs": [],
                    "outputs": {"counts": "${raw.counts}"},
                    "output_artifacts": {
                        "counts": {
                            "type": "bulk.gene-counts",
                            "path": None,
                            "reference": "${raw.counts}",
                            "explicit": override is not None,
                        }
                    },
                    "postflight": postflight,
                    "command_state": "resolved",
                },
                {
                    "id": "downstream",
                    "command": [
                        sys.executable,
                        "-c",
                        "from pathlib import Path; import sys; "
                        "source=Path(sys.argv[1].split('=', 1)[1]); "
                        "Path(sys.argv[2]).write_text(source.read_text())",
                        "input=${raw.counts}",
                        str(downstream_output),
                    ],
                    "cwd": str(workspace),
                    "required_outputs": [str(downstream_output)],
                    "inputs": {"counts": "${raw.counts}"},
                    "outputs": {"copied": str(downstream_output)},
                    "output_artifacts": {
                        "copied": {
                            "type": "test.copy",
                            "path": str(downstream_output),
                            "explicit": True,
                        }
                    },
                    "postflight": {},
                    "command_state": "deferred",
                },
            ]
        }
        return plan, downstream_output

    def test_nfcore_postflight_materializes_deferred_downstream_command(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workspace = root / "work"
            workspace.mkdir()
            plan, downstream_output = self.nfcore_handoff_plan(workspace)
            bundle = self.make_bundle(root, plan)
            result = ProjectExecutor(bundle).execute()

            self.assertTrue(result.succeeded)
            counts = (
                workspace
                / "rnaseq-results"
                / "star_salmon"
                / "salmon.merged.gene_counts.tsv"
            ).resolve()
            self.assertEqual(downstream_output.read_text(), counts.read_text())
            artifact = result.stages[0].artifacts["counts"]
            self.assertEqual(artifact["path"], str(counts))
            self.assertEqual(
                artifact["evidence"]["adapter"], "nfcore-rnaseq-3.26.0"
            )
            self.assertEqual(artifact["evidence"]["release"], "3.26.0")
            self.assertEqual(
                bundle.manifest()["artifacts"]["raw.counts"], artifact
            )
            self.assertIn("downstream.copied", bundle.manifest()["artifacts"])
            downstream_command = result.stages[1].attempts[0].command
            self.assertFalse(
                any("${raw.counts}" in item for item in downstream_command)
            )
            self.assertIn(f"input={counts}", downstream_command)

            resumed = ProjectExecutor(bundle).execute(plan, resume=True)
            self.assertEqual(
                [stage.state for stage in resumed.stages],
                [StageState.SKIPPED, StageState.SKIPPED],
            )
            self.assertEqual(len(resumed.stages[0].attempts), 1)

    def test_nfcore_postflight_honors_explicit_override(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workspace = root / "work"
            workspace.mkdir()
            override = "star_salmon/salmon.merged.gene_counts.tsv"
            plan, downstream_output = self.nfcore_handoff_plan(
                workspace, mode="ambiguous", override=override
            )
            bundle = self.make_bundle(root, plan)
            result = ProjectExecutor(bundle).execute()
            self.assertTrue(result.succeeded)
            self.assertTrue(downstream_output.is_file())
            self.assertEqual(
                result.stages[0].artifacts["counts"]["evidence"]["source"],
                "override",
            )

    def test_nfcore_missing_and_ambiguous_artifacts_fail_postflight(self) -> None:
        cases = (
            ("missing", "No 'bulk.gene-counts'"),
            ("ambiguous", "Ambiguous 'bulk.gene-counts'"),
            ("zero", "violates its output policy"),
        )
        for mode, message in cases:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                workspace = root / "work"
                workspace.mkdir()
                plan, downstream_output = self.nfcore_handoff_plan(
                    workspace, mode=mode
                )
                bundle = self.make_bundle(root, plan)
                result = ProjectExecutor(bundle).execute()
                self.assertEqual(result.status, "failed")
                self.assertEqual(result.stages[0].state, StageState.FAILED)
                self.assertEqual(result.stages[0].attempts[0].exit_code, 0)
                self.assertIn(message, result.stages[0].reason or "")
                self.assertEqual(result.stages[1].state, StageState.SKIPPED)
                self.assertFalse(downstream_output.exists())

    def test_unresolved_deferred_command_fails_without_subprocess(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "must-not-exist"
            plan = {
                "stages": [
                    {
                        "id": "deferred",
                        "command": [
                            sys.executable,
                            "-c",
                            "from pathlib import Path; import sys; "
                            "Path(sys.argv[1]).touch()",
                            str(marker),
                            "${missing.artifact}",
                        ],
                        "command_state": "deferred",
                        "required_outputs": [],
                    }
                ]
            }
            bundle = self.make_bundle(root, plan)
            result = ProjectExecutor(bundle).execute()
            stage = result.stages[0]
            self.assertEqual(stage.state, StageState.FAILED)
            self.assertIsNone(stage.attempts[0].exit_code)
            self.assertIn("unresolved artifact", stage.reason or "")
            self.assertFalse(marker.exists())
            self.assertTrue(bundle.attempt_result_path("deferred", 1).is_file())

    def test_literal_file_content_change_invalidates_resume(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "counts.tsv"
            input_path.write_text("gene\tcount\nA\t1\n", encoding="utf-8")
            output = root / "seen.txt"
            stage = python_stage(
                "consume",
                "from pathlib import Path; import sys; "
                "Path(sys.argv[2]).write_text(Path(sys.argv[1]).read_text())",
                outputs=[str(output)],
                extra_arguments=[str(input_path), str(output)],
                inputs={"counts": str(input_path)},
            )
            stage["input_artifacts"] = {
                "counts": {
                    "type": "bulk.gene-counts",
                    "source": "literal",
                    "value": str(input_path),
                    "references": [],
                    "deferred": False,
                }
            }
            plan = {"stages": [stage]}
            bundle = self.make_bundle(root, plan)
            executor = ProjectExecutor(bundle)
            first = executor.execute()
            initial_hash = first.stages[0].input_fingerprints["counts"]["value"][
                "sha256"
            ]
            self.assertEqual(
                executor.execute(plan, resume=True).stages[0].state,
                StageState.SKIPPED,
            )

            input_path.write_text("gene\tcount\nA\t2\n", encoding="utf-8")
            rerun = executor.execute(plan, resume=True)
            self.assertEqual(rerun.stages[0].state, StageState.COMPLETED)
            self.assertEqual(len(rerun.stages[0].attempts), 2)
            self.assertNotEqual(
                initial_hash,
                rerun.stages[0].input_fingerprints["counts"]["value"]["sha256"],
            )
            attempt = read_json(bundle.attempt_result_path("consume", 2))
            self.assertIn("input_fingerprints", attempt)

    def test_same_path_producer_replacement_invalidates_consumer(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            trigger = root / "trigger.txt"
            trigger.write_text("v1", encoding="utf-8")
            artifact = root / "artifact.tsv"
            copied = root / "copied.tsv"
            producer = python_stage(
                "producer",
                "from pathlib import Path; import sys; "
                "Path(sys.argv[1]).write_text('stable\\n')",
                outputs=[str(artifact)],
                extra_arguments=[str(artifact)],
                inputs={"trigger": str(trigger)},
            )
            producer["input_artifacts"] = {
                "trigger": {
                    "type": "bulk.gene-counts",
                    "value": str(trigger),
                    "references": [],
                }
            }
            producer["output_artifacts"] = {
                "table": {
                    "type": "bulk.gene-counts",
                    "path": str(artifact),
                }
            }
            consumer = python_stage(
                "consumer",
                "from pathlib import Path; import sys; "
                "Path(sys.argv[2]).write_text(Path(sys.argv[1]).read_text())",
                outputs=[str(copied)],
                extra_arguments=["${producer.table}", str(copied)],
                inputs={"table": "${producer.table}"},
            )
            consumer["command_state"] = "deferred"
            consumer["input_artifacts"] = {
                "table": {
                    "type": "bulk.gene-counts",
                    "value": "${producer.table}",
                    "references": [{"stage": "producer", "artifact": "table"}],
                    "deferred": True,
                }
            }
            plan = {"stages": [producer, consumer]}
            bundle = self.make_bundle(root, plan)
            executor = ProjectExecutor(bundle)
            first = executor.execute()
            first_consumer_record = read_json(
                bundle.attempt_result_path("consumer", 1)
            )
            first_fingerprint = first_consumer_record["input_fingerprints"]["table"]
            self.assertTrue(first.succeeded)

            trigger.write_text("v2", encoding="utf-8")
            second = executor.execute(plan, resume=True)
            self.assertEqual(
                [stage.state for stage in second.stages],
                [StageState.COMPLETED, StageState.COMPLETED],
            )
            self.assertEqual(len(second.stages[1].attempts), 2)
            second_fingerprint = second.stages[1].input_fingerprints["table"]
            self.assertEqual(
                first_fingerprint["value"]["sha256"],
                second_fingerprint["value"]["sha256"],
            )
            self.assertNotEqual(
                first_fingerprint["producers"]["producer.table"]["producer"],
                second_fingerprint["producers"]["producer.table"]["producer"],
            )

    def test_known_input_contracts_fail_preflight_without_subprocess(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = (
                ("missing", root / "missing.tsv", "bulk.gene-counts", "does not exist"),
                (
                    "wrong-kind",
                    root / "directory",
                    "bulk.gene-counts",
                    "must be a file",
                ),
                (
                    "bad-sheet",
                    root / "bad.csv",
                    "bulk.rnaseq-samplesheet",
                    "missing columns",
                ),
            )
            (root / "directory").mkdir()
            (root / "bad.csv").write_text(
                "sample,fastq_1\na,a.fastq.gz\n", encoding="utf-8"
            )
            for index, (name, path, artifact_type, message) in enumerate(cases):
                with self.subTest(name=name):
                    marker = root / f"marker-{index}"
                    plan = {
                        "stages": [
                            {
                                "id": name,
                                "command": [
                                    sys.executable,
                                    "-c",
                                    "from pathlib import Path; import sys; "
                                    "Path(sys.argv[1]).touch()",
                                    str(marker),
                                ],
                                "inputs": {"input": str(path)},
                                "input_artifacts": {
                                    "input": {
                                        "type": artifact_type,
                                        "value": str(path),
                                        "references": [],
                                    }
                                },
                                "required_outputs": [],
                            }
                        ]
                    }
                    bundle = self.make_bundle(root, plan, run_id=f"preflight-{index}")
                    stage = ProjectExecutor(bundle).execute().stages[0]
                    self.assertEqual(stage.state, StageState.FAILED)
                    self.assertIsNone(stage.attempts[0].exit_code)
                    self.assertIn(message, stage.reason or "")
                    self.assertFalse(marker.exists())

    def test_stage_outdir_is_created_immediately_before_subprocess(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outdir = root / "fresh" / "nested"
            output = outdir / "ready.txt"
            plan = {
                "stages": [
                    {
                        "id": "outdir",
                        "outdir": str(outdir),
                        "command": [
                            sys.executable,
                            "-c",
                            "from pathlib import Path; import sys; "
                            "root=Path(sys.argv[1]); assert root.is_dir(); "
                            "(root/'ready.txt').write_text('ok')",
                            str(outdir),
                        ],
                        "required_outputs": [str(output)],
                    }
                ]
            }
            bundle = self.make_bundle(root, plan)
            self.assertFalse(outdir.exists())
            result = ProjectExecutor(bundle).execute()
            self.assertTrue(result.succeeded)
            self.assertEqual(output.read_text(), "ok")


if __name__ == "__main__":
    unittest.main()
