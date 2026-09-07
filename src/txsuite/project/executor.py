from __future__ import annotations

import os
import re
import subprocess
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

from txsuite.bulk import validate_samplesheet as validate_bulk_samplesheet
from txsuite.project.adapters.nfcore import resolve_nfcore_artifacts
from txsuite.project.provenance import (
    RunBundle,
    atomic_write_json,
    command_plan_dict,
    fingerprint_path,
    hash_payload,
    read_json,
    redact_arguments,
    utc_now,
    validate_stage_id,
)
from txsuite.project.results import (
    AttemptResult,
    RunResult,
    StageResult,
    StageState,
    decide_resume,
    verify_required_outputs,
)
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import validate_samplesheet as validate_single_cell_samplesheet


_ARTIFACT_REFERENCE = re.compile(
    r"\$\{(?P<stage>[A-Za-z0-9][A-Za-z0-9_-]*)\."
    r"(?P<artifact>[A-Za-z0-9][A-Za-z0-9_-]*)\}"
)
_INPUT_KINDS = {
    "bulk.rnaseq-samplesheet": "file",
    "single-cell.scrnaseq-samplesheet": "file",
    "bulk.gene-counts": "file",
    "sample.metadata": "file",
    "gene-sets.gmt": "file",
    "bulk.differential-expression-table": "file",
    "single-cell.h5ad": "file",
    "single-cell.matrix": "any",
}


class StageExecutionError(RuntimeError):
    def __init__(self, stage: StageResult, run_dir: Path):
        self.stage = stage
        self.run_dir = run_dir
        attempt = stage.attempts[-1] if stage.attempts else None
        detail = attempt.error if attempt is not None else stage.reason
        message = f"Stage {stage.stage_id!r} failed: {detail or 'unknown error'}"
        super().__init__(message)


def _stage_id(stage: Mapping[str, Any]) -> str:
    if "id" not in stage:
        raise ValueError("Every stage must have an 'id'")
    return validate_stage_id(str(stage["id"]))


def select_stages(
    stages: Sequence[Mapping[str, Any]],
    *,
    selected_stage_ids: Sequence[str] | None = None,
    from_stage: str | None = None,
    to_stage: str | None = None,
) -> list[Mapping[str, Any]]:
    """Select stages deterministically, always retaining command-plan order."""
    ids = [_stage_id(stage) for stage in stages]
    if len(ids) != len(set(ids)):
        raise ValueError("Stage IDs must be unique")
    if selected_stage_ids is not None and (
        from_stage is not None or to_stage is not None
    ):
        raise ValueError(
            "selected_stage_ids cannot be combined with from_stage/to_stage"
        )
    if selected_stage_ids is not None:
        selected = {str(item) for item in selected_stage_ids}
        unknown = selected.difference(ids)
        if unknown:
            raise ValueError(f"Unknown selected stage(s): {', '.join(sorted(unknown))}")
        return [stage for stage, stage_id in zip(stages, ids) if stage_id in selected]
    if not ids:
        if from_stage is not None:
            raise ValueError(f"Unknown from_stage: {from_stage}")
        if to_stage is not None:
            raise ValueError(f"Unknown to_stage: {to_stage}")
        return []
    start = ids.index(from_stage) if from_stage is not None and from_stage in ids else 0
    end = (
        ids.index(to_stage)
        if to_stage is not None and to_stage in ids
        else len(ids) - 1
    )
    if from_stage is not None and from_stage not in ids:
        raise ValueError(f"Unknown from_stage: {from_stage}")
    if to_stage is not None and to_stage not in ids:
        raise ValueError(f"Unknown to_stage: {to_stage}")
    if start > end:
        raise ValueError("from_stage must not follow to_stage")
    return list(stages[start : end + 1])


class ProjectExecutor:
    """Execute a planner object or command-plan mapping with durable provenance."""

    def __init__(self, bundle: RunBundle):
        self.bundle = bundle

    def _load_stage_result(self, stage_id: str) -> StageResult | None:
        path = self.bundle.stage_result_path(stage_id)
        if not path.exists():
            return None
        return StageResult.from_dict(read_json(path))

    def _write_stage_result(self, result: StageResult) -> None:
        atomic_write_json(
            self.bundle.stage_result_path(result.stage_id), result.to_dict()
        )
        self.bundle.set_stage_state(result.stage_id, result.state.value)

    @staticmethod
    def _stage_cwd(stage: Mapping[str, Any], default: Path) -> Path:
        cwd_value = stage.get("cwd")
        return Path(str(cwd_value)) if cwd_value is not None else default

    @staticmethod
    def _required_outputs(stage: Mapping[str, Any]) -> tuple[Any, ...]:
        outputs = stage.get("required_outputs", ())
        if isinstance(outputs, (str, Path)) or not isinstance(outputs, Sequence):
            raise ValueError("required_outputs must be a sequence")
        return tuple(outputs)

    @staticmethod
    def _resolved_command(stage: Mapping[str, Any]) -> list[str]:
        command = stage.get("command")
        if isinstance(command, (str, bytes)) or not isinstance(command, Sequence):
            raise ValueError("stage command must be a sequence of arguments")
        resolved = [str(argument) for argument in command]
        if not resolved:
            raise ValueError("stage command must not be empty")
        return resolved

    def _hashes(
        self,
        stage: Mapping[str, Any],
        *,
        input_hash: str | None,
        config_hash: str,
        input_fingerprints: Mapping[str, Any] | None = None,
    ) -> tuple[str, str]:
        declared_input_hash = hash_payload(
            {
                "caller": input_hash,
                "declared": stage.get("inputs", {}),
                "fingerprints": input_fingerprints or {},
            }
        )
        execution_material = {
            "command": stage.get("command"),
            "cwd": stage.get("cwd"),
            "env": stage.get("env", {}),
            "params": stage.get("params", {}),
            "inputs": stage.get("inputs", {}),
            "required_outputs": stage.get("required_outputs", []),
            "input_hash": declared_input_hash,
            "config_hash": config_hash,
        }
        return declared_input_hash, hash_payload(execution_material)

    def _materialize_stage(
        self,
        stage: Mapping[str, Any],
        artifact_map: Mapping[str, Any],
    ) -> dict[str, Any]:
        command = self._resolved_command(stage)
        unresolved: set[str] = set()

        def replace(match: re.Match[str]) -> str:
            key = f"{match.group('stage')}.{match.group('artifact')}"
            artifact = artifact_map.get(key)
            if artifact is None:
                unresolved.add(key)
                return match.group(0)
            if isinstance(artifact, Mapping):
                path = artifact.get("path")
            else:
                path = artifact
            if path is None:
                unresolved.add(key)
                return match.group(0)
            return str(path)

        def materialize(value: Any) -> Any:
            if isinstance(value, str):
                return _ARTIFACT_REFERENCE.sub(replace, value)
            if isinstance(value, Mapping):
                return {str(key): materialize(item) for key, item in value.items()}
            if isinstance(value, (list, tuple)):
                return [materialize(item) for item in value]
            return value

        materialized = [materialize(item) for item in command]
        materialized_inputs = materialize(stage.get("inputs", {}))
        materialized_input_artifacts = materialize(
            stage.get("input_artifacts", {})
        )
        if stage.get("command_state") == "deferred":
            for item in materialized:
                if "${" in item:
                    unresolved.add(item)
        if unresolved:
            names = ", ".join(sorted(unresolved))
            raise ValueError(f"unresolved artifact reference(s): {names}")
        effective = dict(stage)
        effective["command"] = materialized
        effective["inputs"] = materialized_inputs
        if materialized_input_artifacts:
            effective["input_artifacts"] = materialized_input_artifacts
        if effective.get("command_state") == "deferred":
            effective["command_state"] = "resolved"
        return effective

    @staticmethod
    def _path_from_input(value: Any, cwd: Path) -> Path:
        if not isinstance(value, (str, Path)) or not str(value):
            raise ValueError("materialized path input must be a non-empty string")
        path = Path(str(value)).expanduser()
        return path if path.is_absolute() else cwd / path

    def _fingerprint_input_value(self, value: Any, cwd: Path) -> Any:
        if isinstance(value, Mapping):
            return {
                str(key): self._fingerprint_input_value(item, cwd)
                for key, item in value.items()
            }
        if isinstance(value, (list, tuple)):
            return [self._fingerprint_input_value(item, cwd) for item in value]
        if isinstance(value, Path):
            return fingerprint_path(self._path_from_input(value, cwd))
        if isinstance(value, str):
            candidate = self._path_from_input(value, cwd)
            if candidate.exists():
                return fingerprint_path(candidate)
        return {"kind": "literal", "value": value}

    def _input_fingerprints(
        self,
        stage: Mapping[str, Any],
        artifact_map: Mapping[str, Any],
        *,
        cwd: Path,
    ) -> dict[str, Any]:
        inputs = stage.get("inputs", {})
        if not isinstance(inputs, Mapping):
            raise ValueError("stage inputs must be a mapping")
        specifications = stage.get("input_artifacts", {})
        if specifications and not isinstance(specifications, Mapping):
            raise ValueError("stage input_artifacts must be a mapping")
        fingerprints: dict[str, Any] = {}
        for name, value in inputs.items():
            specification = (
                specifications.get(name, {})
                if isinstance(specifications, Mapping)
                else {}
            )
            producers: dict[str, Any] = {}
            if isinstance(specification, Mapping):
                for reference in specification.get("references", []):
                    if not isinstance(reference, Mapping):
                        continue
                    key = f"{reference.get('stage')}.{reference.get('artifact')}"
                    if key in artifact_map:
                        producers[key] = artifact_map[key]
            artifact_type = (
                specification.get("type")
                if isinstance(specification, Mapping)
                else None
            )
            if artifact_type in _INPUT_KINDS:
                value_fingerprint = fingerprint_path(
                    self._path_from_input(value, cwd)
                )
            else:
                value_fingerprint = self._fingerprint_input_value(value, cwd)
            fingerprints[str(name)] = {
                "value": value_fingerprint,
                "artifact_type": artifact_type,
                "producers": producers,
            }
        # File-valued parameters affect execution just as declared inputs do.
        for name in ("metadata", "params_file", "nextflow_config"):
            value = stage.get("params", {}).get(name)
            if value is not None:
                fingerprints[f"params.{name}"] = fingerprint_path(
                    self._path_from_input(value, cwd)
                )
        return fingerprints

    def _validate_materialized_inputs(
        self, stage: Mapping[str, Any], *, cwd: Path
    ) -> None:
        inputs = stage.get("inputs", {})
        specifications = stage.get("input_artifacts", {})
        if not isinstance(inputs, Mapping) or not isinstance(specifications, Mapping):
            return
        for name, specification in specifications.items():
            if not isinstance(specification, Mapping):
                raise ValueError(f"input_artifacts.{name} must be a mapping")
            artifact_type = str(specification.get("type", ""))
            expected_kind = _INPUT_KINDS.get(artifact_type)
            if expected_kind is None:
                continue
            if name not in inputs:
                raise ValueError(f"materialized input {name!r} is missing")
            path = self._path_from_input(inputs[name], cwd)
            if not path.exists():
                raise ValueError(f"input {name!r} does not exist: {path}")
            if expected_kind == "file" and not path.is_file():
                raise ValueError(f"input {name!r} must be a file: {path}")
            if expected_kind == "directory" and not path.is_dir():
                raise ValueError(f"input {name!r} must be a directory: {path}")
            if artifact_type == "bulk.rnaseq-samplesheet":
                validate_bulk_samplesheet(path)
            elif artifact_type == "single-cell.scrnaseq-samplesheet":
                validate_single_cell_samplesheet(path)

    def _resolve_stage_artifacts(
        self, stage: Mapping[str, Any], *, cwd: Path
    ) -> dict[str, dict[str, Any]]:
        stage_id = _stage_id(stage)
        collected: dict[str, dict[str, Any]] = {}
        postflight = stage.get("postflight", {})
        if postflight:
            if not isinstance(postflight, Mapping):
                raise ValueError("stage postflight metadata must be a mapping")
            if postflight.get("adapter") != "nfcore":
                raise ValueError(
                    f"unsupported postflight adapter: {postflight.get('adapter')!r}"
                )
            artifact_types = postflight.get("artifacts", {})
            policies = postflight.get("policies", {})
            overrides = postflight.get("overrides", {})
            if not isinstance(artifact_types, Mapping):
                raise ValueError("postflight artifacts must be a mapping")
            if not isinstance(overrides, Mapping):
                raise ValueError("postflight overrides must be a mapping")
            if not isinstance(policies, Mapping):
                raise ValueError("postflight policies must be a mapping")
            if policies and set(policies) != set(artifact_types):
                raise ValueError(
                    "postflight policies must exactly match declared artifacts"
                )
            results_root = Path(str(postflight["results_root"]))
            if not results_root.is_absolute():
                results_root = cwd / results_root
            for name, artifact_type_value in artifact_types.items():
                artifact_type = str(artifact_type_value)
                override = overrides.get(name)
                resolved = resolve_nfcore_artifacts(
                    results_root,
                    pipeline=str(postflight["pipeline"]),
                    release=str(postflight["release"]),
                    artifact_types=(artifact_type,),
                    overrides=(
                        {artifact_type: str(override)}
                        if override is not None
                        else None
                    ),
                )[artifact_type]
                policy = policies.get(name, {"kind": "any", "non_empty": False})
                if not isinstance(policy, Mapping):
                    raise ValueError(
                        f"postflight policy for {name!r} must be a mapping"
                    )
                verification = verify_required_outputs(
                    [
                        {
                            "path": resolved.path,
                            "kind": policy.get("kind", "any"),
                            "non_empty": policy.get("non_empty", False),
                        }
                    ]
                )
                if not verification.ok:
                    raise ValueError(
                        f"resolved dynamic artifact {name!r} violates its output "
                        f"policy: {', '.join(verification.missing)}"
                    )
                collected[str(name)] = {
                    "stage_id": stage_id,
                    "name": str(name),
                    "artifact_type": artifact_type,
                    "path": str(resolved.path),
                    "evidence": resolved.evidence,
                }

        output_artifacts = stage.get("output_artifacts")
        if isinstance(output_artifacts, Mapping):
            concrete = output_artifacts
        else:
            concrete = {}
        for name, specification in concrete.items():
            if str(name) in collected:
                continue
            if isinstance(specification, Mapping):
                path_value = specification.get("path")
                artifact_type = str(specification.get("type", "unspecified"))
            else:
                path_value = specification
                artifact_type = "unspecified"
            if path_value is None or _ARTIFACT_REFERENCE.search(str(path_value)):
                continue
            path = Path(str(path_value))
            if not path.is_absolute():
                path = cwd / path
            path = path.resolve()
            if not path.exists():
                raise ValueError(
                    f"declared output artifact {stage_id}.{name} does not exist: {path}"
                )
            collected[str(name)] = {
                "stage_id": stage_id,
                "name": str(name),
                "artifact_type": artifact_type,
                "path": str(path),
                "evidence": {"source": "planned"},
            }
        return collected

    @staticmethod
    def _declares_artifacts(stage: Mapping[str, Any]) -> bool:
        postflight = stage.get("postflight")
        if isinstance(postflight, Mapping) and postflight.get("artifacts"):
            return True
        outputs = stage.get("output_artifacts")
        if not isinstance(outputs, Mapping):
            return False
        return any(
            isinstance(specification, Mapping)
            and specification.get("path") is not None
            for specification in outputs.values()
        )

    def _clear_stage_artifacts(
        self, stage_id: str, artifact_map: dict[str, Any]
    ) -> None:
        prefix = f"{stage_id}."
        for key in [key for key in artifact_map if key.startswith(prefix)]:
            del artifact_map[key]
        self.bundle.clear_stage_artifacts(stage_id)

    def _create_stage_outdir(self, stage: Mapping[str, Any]) -> None:
        value = stage.get("outdir")
        if value is None:
            return
        outdir = Path(str(value)).expanduser()
        if not outdir.is_absolute():
            outdir = self._stage_cwd(stage, self.bundle.run_dir) / outdir
        outdir.mkdir(parents=True, exist_ok=True)

    def _persist_attempt(
        self,
        *,
        stage_id: str,
        attempt: AttemptResult,
        input_hash: str,
        config_hash: str,
        execution_hash: str,
        input_fingerprints: Mapping[str, Any],
        required_outputs: tuple[Any, ...],
        artifacts: Mapping[str, Mapping[str, Any]],
    ) -> str:
        path = self.bundle.attempt_result_path(stage_id, attempt.attempt)
        atomic_write_json(
            path,
            {
                "schema_version": 1,
                "run_id": self.bundle.run_id,
                "stage_id": stage_id,
                "input_hash": input_hash,
                "config_hash": config_hash,
                "execution_hash": execution_hash,
                "input_fingerprints": dict(input_fingerprints),
                "required_outputs": list(required_outputs),
                "artifacts": {
                    name: dict(artifact) for name, artifact in artifacts.items()
                },
                "attempt": attempt.to_dict(),
            },
            overwrite=False,
        )
        self.bundle.register_attempt_record(stage_id, path)
        return str(path.relative_to(self.bundle.run_dir))

    def _skipped(
        self,
        stage: Mapping[str, Any],
        *,
        reason: str,
        input_hash: str,
        config_hash: str,
        execution_hash: str,
        input_fingerprints: Mapping[str, Any] | None = None,
        previous: StageResult | None = None,
    ) -> StageResult:
        result = StageResult(
            stage_id=_stage_id(stage),
            state=StageState.SKIPPED,
            input_hash=input_hash,
            config_hash=config_hash,
            execution_hash=execution_hash,
            input_fingerprints=(
                dict(input_fingerprints)
                if input_fingerprints is not None
                else (
                    previous.input_fingerprints if previous is not None else {}
                )
            ),
            required_outputs=self._required_outputs(stage),
            attempts=previous.attempts if previous is not None else (),
            attempt_records=(
                previous.attempt_records if previous is not None else ()
            ),
            artifacts=previous.artifacts if previous is not None else {},
            reason=reason,
            updated_at=utc_now(),
        )
        self._write_stage_result(result)
        return result

    def _execute_stage(
        self,
        stage: Mapping[str, Any],
        *,
        input_hash: str,
        config_hash: str,
        execution_hash: str,
        input_fingerprints: Mapping[str, Any],
        previous: StageResult | None,
    ) -> StageResult:
        stage_id = _stage_id(stage)
        command = self._resolved_command(stage)
        required_outputs = self._required_outputs(stage)
        cwd = self._stage_cwd(stage, self.bundle.run_dir)
        env_value = stage.get("env", {})
        if not isinstance(env_value, Mapping):
            raise ValueError("stage env must be a mapping")
        environment = os.environ.copy()
        environment.update({str(key): str(value) for key, value in env_value.items()})

        old_attempts = previous.attempts if previous is not None else ()
        attempt_number = len(old_attempts) + 1
        stdout_path = self.bundle.log_path(stage_id, attempt_number, "stdout")
        stderr_path = self.bundle.log_path(stage_id, attempt_number, "stderr")
        stdout_relative = str(stdout_path.relative_to(self.bundle.run_dir))
        stderr_relative = str(stderr_path.relative_to(self.bundle.run_dir))
        started_at = utc_now()
        started = time.monotonic()
        running_attempt = AttemptResult(
            attempt=attempt_number,
            state=StageState.RUNNING,
            command=tuple(redact_arguments(command)),
            stdout_log=stdout_relative,
            stderr_log=stderr_relative,
            started_at=started_at,
        )
        running = StageResult(
            stage_id=stage_id,
            state=StageState.RUNNING,
            input_hash=input_hash,
            config_hash=config_hash,
            execution_hash=execution_hash,
            input_fingerprints=dict(input_fingerprints),
            required_outputs=required_outputs,
            attempts=(*old_attempts, running_attempt),
            attempt_records=(
                previous.attempt_records if previous is not None else ()
            ),
            updated_at=started_at,
        )
        self._write_stage_result(running)

        exit_code: int | None = None
        error: str | None = None
        verification = None
        artifacts: dict[str, dict[str, Any]] = {}
        try:
            with (
                stdout_path.open("w", encoding="utf-8") as stdout,
                stderr_path.open("w", encoding="utf-8") as stderr,
            ):
                completed = subprocess.run(
                    command,
                    cwd=cwd,
                    env=environment,
                    stdout=stdout,
                    stderr=stderr,
                    text=True,
                    check=False,
                    shell=False,
                )
            exit_code = completed.returncode
            if exit_code != 0:
                error = f"command exited with status {exit_code}"
            else:
                verification = verify_required_outputs(required_outputs, base_dir=cwd)
                if not verification.ok:
                    error = "required outputs are missing or invalid: " + ", ".join(
                        verification.missing
                    )
                else:
                    try:
                        artifacts = self._resolve_stage_artifacts(stage, cwd=cwd)
                    except (KeyError, TypeError, ValueError, TxSuiteError) as exc:
                        error = f"postflight artifact resolution failed: {exc}"
        except OSError as exc:
            error = f"cannot execute {command[0]!r}: {exc}"

        finished_at = utc_now()
        state = StageState.COMPLETED if error is None else StageState.FAILED
        if state is StageState.COMPLETED:
            for artifact in artifacts.values():
                artifact["producer"] = {
                    "run_id": self.bundle.run_id,
                    "execution_hash": execution_hash,
                    "attempt": attempt_number,
                }
        finished_attempt = AttemptResult(
            attempt=attempt_number,
            state=state,
            command=tuple(redact_arguments(command)),
            stdout_log=stdout_relative,
            stderr_log=stderr_relative,
            started_at=started_at,
            finished_at=finished_at,
            duration_seconds=max(0.0, time.monotonic() - started),
            exit_code=exit_code,
            error=error,
            output_verification=verification,
        )
        attempt_record = self._persist_attempt(
            stage_id=stage_id,
            attempt=finished_attempt,
            input_hash=input_hash,
            config_hash=config_hash,
            execution_hash=execution_hash,
            input_fingerprints=input_fingerprints,
            required_outputs=required_outputs,
            artifacts=artifacts,
        )
        old_attempt_records = (
            previous.attempt_records if previous is not None else ()
        )
        if state is StageState.COMPLETED and artifacts:
            self.bundle.record_artifacts(stage_id, artifacts)
        result = StageResult(
            stage_id=stage_id,
            state=state,
            input_hash=input_hash,
            config_hash=config_hash,
            execution_hash=execution_hash,
            input_fingerprints=dict(input_fingerprints),
            required_outputs=required_outputs,
            attempts=(*old_attempts, finished_attempt),
            attempt_records=(*old_attempt_records, attempt_record),
            artifacts=artifacts,
            reason=error,
            updated_at=finished_at,
        )
        self._write_stage_result(result)
        return result

    def _record_preflight_failure(
        self,
        stage: Mapping[str, Any],
        *,
        error: str,
        input_hash: str,
        config_hash: str,
        execution_hash: str,
        input_fingerprints: Mapping[str, Any] | None = None,
        previous: StageResult | None,
    ) -> StageResult:
        stage_id = _stage_id(stage)
        command = self._resolved_command(stage)
        required_outputs = self._required_outputs(stage)
        old_attempts = previous.attempts if previous is not None else ()
        old_records = previous.attempt_records if previous is not None else ()
        attempt_number = len(old_attempts) + 1
        stdout_path = self.bundle.log_path(stage_id, attempt_number, "stdout")
        stderr_path = self.bundle.log_path(stage_id, attempt_number, "stderr")
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text("", encoding="utf-8")
        timestamp = utc_now()
        attempt = AttemptResult(
            attempt=attempt_number,
            state=StageState.FAILED,
            command=tuple(redact_arguments(command)),
            stdout_log=str(stdout_path.relative_to(self.bundle.run_dir)),
            stderr_log=str(stderr_path.relative_to(self.bundle.run_dir)),
            started_at=timestamp,
            finished_at=timestamp,
            duration_seconds=0.0,
            exit_code=None,
            error=f"preflight failed: {error}",
        )
        record = self._persist_attempt(
            stage_id=stage_id,
            attempt=attempt,
            input_hash=input_hash,
            config_hash=config_hash,
            execution_hash=execution_hash,
            input_fingerprints=input_fingerprints or {},
            required_outputs=required_outputs,
            artifacts={},
        )
        result = StageResult(
            stage_id=stage_id,
            state=StageState.FAILED,
            input_hash=input_hash,
            config_hash=config_hash,
            execution_hash=execution_hash,
            input_fingerprints=dict(input_fingerprints or {}),
            required_outputs=required_outputs,
            attempts=(*old_attempts, attempt),
            attempt_records=(*old_records, record),
            artifacts={},
            reason=attempt.error,
            updated_at=timestamp,
        )
        self._write_stage_result(result)
        return result

    def execute(
        self,
        command_plan: Any | None = None,
        *,
        resume: bool = False,
        input_hashes: Mapping[str, str] | None = None,
        config_hash: str | None = None,
        selected_stage_ids: Sequence[str] | None = None,
        from_stage: str | None = None,
        to_stage: str | None = None,
        fail_fast: bool = True,
        raise_on_failure: bool = False,
    ) -> RunResult:
        plan = (
            command_plan_dict(command_plan)
            if command_plan is not None
            else self.bundle.command_plan()
        )
        stages_value = plan.get("stages", [])
        if not isinstance(stages_value, list):
            raise ValueError("command_plan['stages'] must be a list")
        stages: list[Mapping[str, Any]] = []
        for stage in stages_value:
            if not isinstance(stage, Mapping):
                raise ValueError("command-plan stages must be mappings")
            stages.append(stage)
        chosen = select_stages(
            stages,
            selected_stage_ids=selected_stage_ids,
            from_stage=from_stage,
            to_stage=to_stage,
        )
        chosen_ids = {_stage_id(stage) for stage in chosen}
        resolved_config_hash = config_hash or str(
            self.bundle.manifest().get("hashes", {}).get(
                "resolved_config", hash_payload(self.bundle.resolved_config())
            )
        )
        supplied_input_hashes = input_hashes or {}
        started_at = utc_now()
        self.bundle.update_manifest(status="running", started_at=started_at)
        results: list[StageResult] = []
        failed: StageResult | None = None
        artifact_map: dict[str, Any] = dict(
            self.bundle.manifest().get("artifacts", {})
        )

        for stage in stages:
            stage_id = _stage_id(stage)
            original_input_hash, original_execution_hash = self._hashes(
                stage,
                input_hash=supplied_input_hashes.get(stage_id),
                config_hash=resolved_config_hash,
            )
            previous = self._load_stage_result(stage_id)
            if stage_id not in chosen_ids:
                if previous is not None:
                    results.append(previous)
                else:
                    results.append(
                        self._skipped(
                            stage,
                            reason="not selected",
                            input_hash=original_input_hash,
                            config_hash=resolved_config_hash,
                            execution_hash=original_execution_hash,
                        )
                    )
                continue
            if failed is not None and fail_fast:
                results.append(
                    self._skipped(
                        stage,
                        reason=f"blocked by failed stage {failed.stage_id}",
                        input_hash=original_input_hash,
                        config_hash=resolved_config_hash,
                        execution_hash=original_execution_hash,
                        previous=previous,
                    )
                )
                continue
            try:
                effective_stage = self._materialize_stage(stage, artifact_map)
            except ValueError as exc:
                self._clear_stage_artifacts(stage_id, artifact_map)
                result = self._record_preflight_failure(
                    stage,
                    error=str(exc),
                    input_hash=original_input_hash,
                    config_hash=resolved_config_hash,
                    execution_hash=original_execution_hash,
                    previous=previous,
                )
                results.append(result)
                if failed is None:
                    failed = result
                continue
            input_cwd = self._stage_cwd(effective_stage, self.bundle.run_dir)
            input_fingerprints: dict[str, Any] = {}
            try:
                input_fingerprints = self._input_fingerprints(
                    effective_stage, artifact_map, cwd=input_cwd
                )
                self._validate_materialized_inputs(
                    effective_stage, cwd=input_cwd
                )
            except (OSError, TypeError, ValueError, TxSuiteError) as exc:
                input_fingerprints["_validation_error"] = str(exc)
                failed_input_hash, failed_execution_hash = self._hashes(
                    effective_stage,
                    input_hash=supplied_input_hashes.get(stage_id),
                    config_hash=resolved_config_hash,
                    input_fingerprints=input_fingerprints,
                )
                self._clear_stage_artifacts(stage_id, artifact_map)
                result = self._record_preflight_failure(
                    effective_stage,
                    error=f"input validation failed: {exc}",
                    input_hash=failed_input_hash,
                    config_hash=resolved_config_hash,
                    execution_hash=failed_execution_hash,
                    input_fingerprints=input_fingerprints,
                    previous=previous,
                )
                results.append(result)
                if failed is None:
                    failed = result
                continue
            declared_hash, execution_hash = self._hashes(
                effective_stage,
                input_hash=supplied_input_hashes.get(stage_id),
                config_hash=resolved_config_hash,
                input_fingerprints=input_fingerprints,
            )
            if resume:
                output_base = self._stage_cwd(
                    effective_stage, self.bundle.run_dir
                )
                resume_outputs = list(self._required_outputs(effective_stage))
                if previous is not None:
                    resume_outputs.extend(
                        artifact["path"]
                        for artifact in previous.artifacts.values()
                        if "path" in artifact
                    )
                decision = decide_resume(
                    previous,
                    input_hash=declared_hash,
                    config_hash=resolved_config_hash,
                    execution_hash=execution_hash,
                    required_outputs=resume_outputs,
                    base_dir=output_base,
                )
                has_required_artifact_provenance = (
                    not self._declares_artifacts(effective_stage)
                    or (previous is not None and bool(previous.artifacts))
                )
                if decision.should_skip and has_required_artifact_provenance:
                    results.append(
                        self._skipped(
                            effective_stage,
                            reason=f"resume: {decision.reason}",
                            input_hash=declared_hash,
                            config_hash=resolved_config_hash,
                            execution_hash=execution_hash,
                            input_fingerprints=input_fingerprints,
                            previous=previous,
                        )
                    )
                    continue
            self._clear_stage_artifacts(stage_id, artifact_map)
            try:
                self._create_stage_outdir(effective_stage)
            except OSError as exc:
                result = self._record_preflight_failure(
                    effective_stage,
                    error=f"cannot create stage outdir: {exc}",
                    input_hash=declared_hash,
                    config_hash=resolved_config_hash,
                    execution_hash=execution_hash,
                    input_fingerprints=input_fingerprints,
                    previous=previous,
                )
                results.append(result)
                if failed is None:
                    failed = result
                continue
            result = self._execute_stage(
                effective_stage,
                input_hash=declared_hash,
                config_hash=resolved_config_hash,
                execution_hash=execution_hash,
                input_fingerprints=input_fingerprints,
                previous=previous,
            )
            results.append(result)
            if result.state is StageState.COMPLETED:
                for name, artifact in result.artifacts.items():
                    artifact_map[f"{stage_id}.{name}"] = artifact
            if result.state is StageState.FAILED and failed is None:
                failed = result

        finished_at = utc_now()
        status = "failed" if failed is not None else "completed"
        self.bundle.update_manifest(status=status, finished_at=finished_at)
        run_result = RunResult(
            run_id=self.bundle.run_id,
            status=status,
            stages=tuple(results),
            started_at=started_at,
            finished_at=finished_at,
        )
        if failed is not None and raise_on_failure:
            raise StageExecutionError(failed, self.bundle.run_dir)
        return run_result
