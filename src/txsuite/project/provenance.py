from __future__ import annotations

import hashlib
import json
import os
import re
import secrets
import tempfile
from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = 1
# Keep provenance snapshots consistent with project-schema redaction output.
REDACTED = "***REDACTED***"
_SAFE_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
_STAGE_STATES = frozenset(
    {"pending", "running", "completed", "failed", "skipped"}
)
_SENSITIVE_NAMES = {
    "access-key",
    "api-key",
    "apikey",
    "auth",
    "authorization",
    "credential",
    "credentials",
    "passphrase",
    "passwd",
    "password",
    "private-key",
    "secret",
    "token",
}


class ProvenanceError(RuntimeError):
    """Raised when a provenance bundle cannot be created or updated safely."""


def utc_now() -> str:
    return datetime.now(UTC).isoformat()


def generate_run_id(prefix: str = "run") -> str:
    """Return a sortable run identifier with enough entropy for concurrent callers."""
    if not _SAFE_NAME.fullmatch(prefix):
        raise ValueError(f"Unsafe run ID prefix: {prefix!r}")
    timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S.%fZ")
    return f"{prefix}-{timestamp}-{secrets.token_hex(4)}"


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if hasattr(value, "to_dict"):
        return value.to_dict()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def canonical_json(payload: Any) -> str:
    return json.dumps(
        payload,
        default=_json_default,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def hash_payload(payload: Any) -> str:
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint_path(path: Path) -> dict[str, Any]:
    """Return a deterministic, content-addressed fingerprint for a filesystem path."""
    path = Path(path)
    absolute = path.resolve(strict=False)
    if not path.exists():
        return {"path": str(absolute), "kind": "missing"}
    if path.is_file():
        return {
            "path": str(absolute),
            "kind": "file",
            "size": path.stat().st_size,
            "sha256": hash_file(path),
        }
    if path.is_dir():
        entries: list[dict[str, Any]] = []
        for child in sorted(
            path.rglob("*"), key=lambda item: item.relative_to(path).as_posix()
        ):
            relative = child.relative_to(path).as_posix()
            if child.is_symlink():
                entry = {
                    "path": relative,
                    "kind": "symlink",
                    "target": os.readlink(child),
                }
            elif child.is_file():
                entry = {
                    "path": relative,
                    "kind": "file",
                    "size": child.stat().st_size,
                    "sha256": hash_file(child),
                }
            elif child.is_dir():
                entry = {"path": relative, "kind": "directory"}
            else:
                entry = {"path": relative, "kind": "other"}
            entries.append(entry)
        return {
            "path": str(absolute),
            "kind": "directory",
            "sha256": hash_payload(entries),
            "entries": entries,
        }
    return {"path": str(absolute), "kind": "other"}


def atomic_write_json(
    path: Path, payload: Any, *, overwrite: bool = True
) -> None:
    """Durably replace a JSON file without exposing a partially written document."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not overwrite and path.exists():
        raise FileExistsError(path)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            json.dump(
                payload,
                handle,
                default=_json_default,
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if overwrite:
            os.replace(temporary_name, path)
            temporary_name = None
        else:
            # Linking is an atomic no-clobber publish on the same filesystem.
            try:
                os.link(temporary_name, path)
            except FileExistsError:
                raise FileExistsError(path) from None
            Path(temporary_name).unlink()
            temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ProvenanceError(f"Expected a JSON object in {path}")
    return payload


def _sensitive_name(value: str) -> bool:
    normalized = value.lower().lstrip("-").replace("_", "-")
    return normalized in _SENSITIVE_NAMES or any(
        normalized.endswith(f"-{name}") for name in _SENSITIVE_NAMES
    )


def _redact_url_userinfo(value: str) -> str:
    return re.sub(r"(://)[^/@\s]+(?::[^/@\s]*)?@", rf"\1{REDACTED}@", value)


def redact_arguments(command: Sequence[str]) -> list[str]:
    """Redact common secret-bearing CLI forms while retaining a useful command trace."""
    redacted: list[str] = []
    redact_next = False
    for raw in command:
        argument = str(raw)
        if redact_next:
            redacted.append(REDACTED)
            redact_next = False
            continue
        if "=" in argument:
            name, value = argument.split("=", 1)
            if _sensitive_name(name):
                redacted.append(f"{name}={REDACTED}")
                continue
        if argument.startswith("-") and _sensitive_name(argument):
            redacted.append(argument)
            redact_next = True
            continue
        redacted.append(_redact_url_userinfo(argument))
    return redacted


def _redact_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return redact_mapping(value)
    if isinstance(value, (list, tuple)):
        return [_redact_value(item) for item in value]
    return value


def redact_mapping(payload: Mapping[str, Any]) -> dict[str, Any]:
    redacted: dict[str, Any] = {}
    for key, value in payload.items():
        if _sensitive_name(str(key)):
            redacted[str(key)] = REDACTED
        else:
            redacted[str(key)] = _redact_value(value)
    return redacted


def redact_command_plan(command_plan: Mapping[str, Any]) -> dict[str, Any]:
    plan = redact_mapping(command_plan)
    stages = plan.get("stages", [])
    if isinstance(stages, list):
        for stage in stages:
            command = stage.get("command") if isinstance(stage, dict) else None
            if isinstance(command, Sequence) and not isinstance(command, (str, bytes)):
                stage["command"] = redact_arguments(command)
    return plan


def validate_stage_id(stage_id: str) -> str:
    if not _SAFE_NAME.fullmatch(stage_id):
        raise ValueError(f"Unsafe stage ID: {stage_id!r}")
    return stage_id


def command_plan_dict(command_plan: Any) -> dict[str, Any]:
    """Normalize a planner object or an already serialized command plan."""
    value = command_plan
    if not isinstance(value, Mapping):
        to_dict = getattr(value, "to_dict", None)
        if not callable(to_dict):
            raise TypeError("command_plan must be a mapping or provide to_dict()")
        value = to_dict()
    if not isinstance(value, Mapping):
        raise TypeError("command_plan.to_dict() must return a mapping")
    return dict(value)


def _command_plan_stage_ids(command_plan: Mapping[str, Any]) -> list[str]:
    stages = command_plan.get("stages", [])
    if not isinstance(stages, list):
        raise ValueError("command_plan['stages'] must be a list")
    stage_ids: list[str] = []
    for stage in stages:
        if not isinstance(stage, Mapping) or "id" not in stage:
            raise ValueError("Every command-plan stage must have an 'id'")
        command = stage.get("command")
        if (
            isinstance(command, (str, bytes))
            or not isinstance(command, Sequence)
            or not command
        ):
            raise ValueError("Every command-plan stage needs a non-empty command list")
        outputs = stage.get("required_outputs", [])
        if isinstance(outputs, (str, bytes)) or not isinstance(outputs, Sequence):
            raise ValueError("Stage required_outputs must be a list")
        if not isinstance(stage.get("env", {}), Mapping):
            raise ValueError("Stage env must be a mapping")
        stage_ids.append(validate_stage_id(str(stage["id"])))
    if len(stage_ids) != len(set(stage_ids)):
        raise ValueError("Command-plan stage IDs must be unique")
    return stage_ids


@dataclass(frozen=True)
class RunBundle:
    """Filesystem-backed provenance for a single project run."""

    run_dir: Path
    run_id: str
    _command_plan: Mapping[str, Any] | None = field(
        default=None, repr=False, compare=False
    )

    @classmethod
    def create(
        cls,
        runs_dir: Path,
        *,
        workflow: Mapping[str, Any],
        resolved_config: Mapping[str, Any],
        command_plan: Any,
        run_id: str | None = None,
    ) -> "RunBundle":
        # Validate all obvious caller-controlled names before allocating a run
        # directory. A partial directory then means an interrupted write, not a
        # predictable input error.
        plan_snapshot = command_plan_dict(command_plan)
        stage_ids = _command_plan_stage_ids(plan_snapshot)
        identifier = run_id or generate_run_id()
        if not _SAFE_NAME.fullmatch(identifier):
            raise ValueError(f"Unsafe run ID: {identifier!r}")

        workflow_snapshot = dict(workflow)
        config_snapshot = dict(resolved_config)
        workflow_hash = hash_payload(workflow_snapshot)
        config_hash = hash_payload(config_snapshot)
        plan_hash = hash_payload(plan_snapshot)
        redacted_workflow = redact_mapping(workflow_snapshot)
        redacted_config = redact_mapping(config_snapshot)
        redacted_plan = redact_command_plan(plan_snapshot)

        runs_dir = Path(runs_dir)
        runs_dir.mkdir(parents=True, exist_ok=True)
        run_dir = runs_dir / identifier
        try:
            run_dir.mkdir()
            (run_dir / "logs").mkdir()
            (run_dir / "stage-results").mkdir()
        except FileExistsError as exc:
            raise ProvenanceError(f"Run bundle already exists: {run_dir}") from exc

        created_at = utc_now()
        try:
            atomic_write_json(
                run_dir / "workflow.json",
                redacted_workflow,
                overwrite=False,
            )
            atomic_write_json(
                run_dir / "resolved-config.json",
                redacted_config,
                overwrite=False,
            )
            atomic_write_json(
                run_dir / "command-plan.json",
                redacted_plan,
                overwrite=False,
            )
            atomic_write_json(
                run_dir / "run-manifest.json",
                {
                    "schema_version": SCHEMA_VERSION,
                    "run_id": identifier,
                    "status": "pending",
                    "created_at": created_at,
                    "updated_at": created_at,
                    "files": {
                        "workflow": "workflow.json",
                        "resolved_config": "resolved-config.json",
                        "command_plan": "command-plan.json",
                        "logs": "logs",
                        "stage_results": "stage-results",
                    },
                    "hashes": {
                        "workflow": workflow_hash,
                        "resolved_config": config_hash,
                        "command_plan": plan_hash,
                    },
                    "stage_states": {
                        stage_id: "pending" for stage_id in stage_ids
                    },
                    "stage_attempt_records": {
                        stage_id: [] for stage_id in stage_ids
                    },
                    "artifacts": {},
                },
                overwrite=False,
            )
        except Exception:
            # The exclusive run directory remains as evidence of an interrupted create;
            # callers must choose a new run ID rather than mutating a partial bundle.
            raise
        return cls(run_dir=run_dir, run_id=identifier, _command_plan=plan_snapshot)

    @classmethod
    def load(cls, run_dir: Path) -> "RunBundle":
        run_dir = Path(run_dir)
        manifest = read_json(run_dir / "run-manifest.json")
        return cls(run_dir=run_dir, run_id=str(manifest["run_id"]))

    @property
    def manifest_path(self) -> Path:
        return self.run_dir / "run-manifest.json"

    def manifest(self) -> dict[str, Any]:
        return read_json(self.manifest_path)

    def command_plan(self) -> dict[str, Any]:
        if self._command_plan is not None:
            return dict(self._command_plan)
        return read_json(self.run_dir / "command-plan.json")

    def resolved_config(self) -> dict[str, Any]:
        return read_json(self.run_dir / "resolved-config.json")

    def update_manifest(self, **changes: Any) -> dict[str, Any]:
        manifest = self.manifest()
        manifest.update(changes)
        manifest["updated_at"] = utc_now()
        atomic_write_json(self.manifest_path, manifest)
        return manifest

    def set_stage_state(self, stage_id: str, state: str) -> None:
        validate_stage_id(stage_id)
        state_value = str(getattr(state, "value", state))
        if state_value not in _STAGE_STATES:
            raise ValueError(f"Unknown stage state: {state_value!r}")
        manifest = self.manifest()
        states = manifest.setdefault("stage_states", {})
        states[stage_id] = state_value
        manifest["updated_at"] = utc_now()
        atomic_write_json(self.manifest_path, manifest)

    def stage_result_path(self, stage_id: str) -> Path:
        return self.run_dir / "stage-results" / f"{validate_stage_id(stage_id)}.json"

    def attempt_result_path(self, stage_id: str, attempt: int) -> Path:
        if attempt < 1:
            raise ValueError("attempt must be positive")
        name = (
            f"{validate_stage_id(stage_id)}.attempt-{attempt}.{self.run_id}.json"
        )
        return self.run_dir / "stage-results" / name

    def register_attempt_record(self, stage_id: str, record_path: Path) -> None:
        validate_stage_id(stage_id)
        relative = str(Path(record_path).relative_to(self.run_dir))
        manifest = self.manifest()
        records = manifest.setdefault("stage_attempt_records", {}).setdefault(
            stage_id, []
        )
        if relative not in records:
            records.append(relative)
        manifest["updated_at"] = utc_now()
        atomic_write_json(self.manifest_path, manifest)

    def record_artifacts(
        self, stage_id: str, artifacts: Mapping[str, Mapping[str, Any]]
    ) -> None:
        validate_stage_id(stage_id)
        manifest = self.manifest()
        recorded = manifest.setdefault("artifacts", {})
        for name, artifact in artifacts.items():
            recorded[f"{stage_id}.{name}"] = dict(artifact)
        manifest["updated_at"] = utc_now()
        atomic_write_json(self.manifest_path, manifest)

    def clear_stage_artifacts(self, stage_id: str) -> None:
        validate_stage_id(stage_id)
        manifest = self.manifest()
        recorded = manifest.setdefault("artifacts", {})
        prefix = f"{stage_id}."
        for key in [key for key in recorded if key.startswith(prefix)]:
            del recorded[key]
        manifest["updated_at"] = utc_now()
        atomic_write_json(self.manifest_path, manifest)

    def log_path(self, stage_id: str, attempt: int, stream: str) -> Path:
        if stream not in {"stdout", "stderr"}:
            raise ValueError("stream must be 'stdout' or 'stderr'")
        if attempt < 1:
            raise ValueError("attempt must be positive")
        name = f"{validate_stage_id(stage_id)}.attempt-{attempt}.{stream}.log"
        return self.run_dir / "logs" / name
