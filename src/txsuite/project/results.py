from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence


class StageState(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


@dataclass(frozen=True)
class OutputVerification:
    ok: bool
    checks: tuple[dict[str, Any], ...] = ()

    @property
    def missing(self) -> tuple[str, ...]:
        return tuple(str(check["path"]) for check in self.checks if not check["valid"])

    def to_dict(self) -> dict[str, Any]:
        return {"ok": self.ok, "checks": [dict(check) for check in self.checks]}


def verify_required_outputs(
    outputs: Sequence[str | Path | Mapping[str, Any]],
    *,
    base_dir: Path | None = None,
) -> OutputVerification:
    """Verify required paths, optionally enforcing type and non-empty constraints."""
    checks: list[dict[str, Any]] = []
    root = Path(base_dir) if base_dir is not None else None
    for output in outputs:
        if isinstance(output, Mapping):
            if "path" not in output:
                raise ValueError("Required-output mappings need a 'path'")
            raw_path = Path(str(output["path"]))
            expected_kind = str(output.get("kind", "any"))
            non_empty = bool(output.get("non_empty", False))
        else:
            raw_path = Path(output)
            expected_kind = "any"
            non_empty = False
        path = raw_path if raw_path.is_absolute() or root is None else root / raw_path
        exists = path.exists()
        kind_ok = (
            exists
            and (
                expected_kind == "any"
                or (expected_kind == "file" and path.is_file())
                or (expected_kind in {"dir", "directory"} and path.is_dir())
            )
        )
        if expected_kind not in {"any", "file", "dir", "directory"}:
            raise ValueError(f"Unknown required-output kind: {expected_kind!r}")
        size_ok = True
        if kind_ok and non_empty:
            if path.is_file():
                size_ok = path.stat().st_size > 0
            elif path.is_dir():
                size_ok = next(path.iterdir(), None) is not None
        valid = bool(kind_ok and size_ok)
        reason = None
        if not exists:
            reason = "missing"
        elif not kind_ok:
            reason = f"expected {expected_kind}"
        elif not size_ok:
            reason = "empty"
        checks.append(
            {
                "path": str(path),
                "kind": expected_kind,
                "non_empty": non_empty,
                "valid": valid,
                "reason": reason,
            }
        )
    return OutputVerification(
        ok=all(check["valid"] for check in checks), checks=tuple(checks)
    )


@dataclass(frozen=True)
class AttemptResult:
    attempt: int
    state: StageState
    command: tuple[str, ...]
    stdout_log: str
    stderr_log: str
    started_at: str
    finished_at: str | None = None
    duration_seconds: float | None = None
    exit_code: int | None = None
    error: str | None = None
    output_verification: OutputVerification | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "attempt": self.attempt,
            "state": self.state.value,
            "command": list(self.command),
            "stdout_log": self.stdout_log,
            "stderr_log": self.stderr_log,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "duration_seconds": self.duration_seconds,
            "exit_code": self.exit_code,
            "error": self.error,
            "output_verification": (
                self.output_verification.to_dict()
                if self.output_verification is not None
                else None
            ),
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "AttemptResult":
        verification = payload.get("output_verification")
        return cls(
            attempt=int(payload["attempt"]),
            state=StageState(str(payload["state"])),
            command=tuple(str(item) for item in payload.get("command", [])),
            stdout_log=str(payload.get("stdout_log", "")),
            stderr_log=str(payload.get("stderr_log", "")),
            started_at=str(payload["started_at"]),
            finished_at=payload.get("finished_at"),
            duration_seconds=payload.get("duration_seconds"),
            exit_code=payload.get("exit_code"),
            error=payload.get("error"),
            output_verification=(
                OutputVerification(
                    ok=bool(verification["ok"]),
                    checks=tuple(dict(item) for item in verification.get("checks", [])),
                )
                if isinstance(verification, Mapping)
                else None
            ),
        )


@dataclass(frozen=True)
class StageResult:
    stage_id: str
    state: StageState
    input_hash: str
    config_hash: str
    execution_hash: str
    required_outputs: tuple[Any, ...] = ()
    attempts: tuple[AttemptResult, ...] = ()
    attempt_records: tuple[str, ...] = ()
    artifacts: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)
    reason: str | None = None
    updated_at: str | None = None
    input_fingerprints: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": 1,
            "stage_id": self.stage_id,
            "state": self.state.value,
            "input_hash": self.input_hash,
            "config_hash": self.config_hash,
            "execution_hash": self.execution_hash,
            "input_fingerprints": dict(self.input_fingerprints),
            "required_outputs": list(self.required_outputs),
            "attempts": [attempt.to_dict() for attempt in self.attempts],
            "attempt_records": list(self.attempt_records),
            "artifacts": {
                str(name): dict(artifact)
                for name, artifact in self.artifacts.items()
            },
            "reason": self.reason,
            "updated_at": self.updated_at,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "StageResult":
        return cls(
            stage_id=str(payload["stage_id"]),
            state=StageState(str(payload["state"])),
            input_hash=str(payload.get("input_hash", "")),
            config_hash=str(payload.get("config_hash", "")),
            execution_hash=str(payload.get("execution_hash", "")),
            input_fingerprints=dict(payload.get("input_fingerprints", {})),
            required_outputs=tuple(payload.get("required_outputs", [])),
            attempts=tuple(
                AttemptResult.from_dict(item) for item in payload.get("attempts", [])
            ),
            attempt_records=tuple(
                str(item) for item in payload.get("attempt_records", [])
            ),
            artifacts={
                str(name): dict(artifact)
                for name, artifact in payload.get("artifacts", {}).items()
            },
            reason=payload.get("reason"),
            updated_at=payload.get("updated_at"),
        )


@dataclass(frozen=True)
class ResumeDecision:
    should_run: bool
    reason: str
    verification: OutputVerification | None = None

    @property
    def should_skip(self) -> bool:
        return not self.should_run


def decide_resume(
    previous: StageResult | Mapping[str, Any] | None,
    *,
    input_hash: str,
    config_hash: str,
    execution_hash: str,
    required_outputs: Sequence[str | Path | Mapping[str, Any]] = (),
    base_dir: Path | None = None,
) -> ResumeDecision:
    if previous is None:
        return ResumeDecision(True, "no previous stage result")
    result = (
        previous
        if isinstance(previous, StageResult)
        else StageResult.from_dict(previous)
    )
    resume_skip = result.state is StageState.SKIPPED and bool(
        result.reason and result.reason.startswith("resume:")
    )
    if result.state is not StageState.COMPLETED and not resume_skip:
        return ResumeDecision(True, f"previous state is {result.state.value}")
    if result.input_hash != input_hash:
        return ResumeDecision(True, "declared inputs changed")
    if result.config_hash != config_hash:
        return ResumeDecision(True, "resolved config changed")
    if result.execution_hash != execution_hash:
        return ResumeDecision(True, "resolved stage command or parameters changed")
    verification = verify_required_outputs(required_outputs, base_dir=base_dir)
    if not verification.ok:
        return ResumeDecision(
            True, "required outputs are missing or invalid", verification
        )
    return ResumeDecision(
        False, "input, config, command, and output checks match", verification
    )


@dataclass(frozen=True)
class RunResult:
    run_id: str
    status: str
    stages: tuple[StageResult, ...] = field(default_factory=tuple)
    started_at: str | None = None
    finished_at: str | None = None

    @property
    def succeeded(self) -> bool:
        return self.status == "completed"

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "status": self.status,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "stages": [stage.to_dict() for stage in self.stages],
        }
