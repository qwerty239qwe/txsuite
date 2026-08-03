"""Load and resolve strict TxSuite project workflow TOML files."""

from __future__ import annotations

import math
import tomllib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from txsuite.project.model import (
    ArtifactReference,
    ExecutionSettings,
    ProjectDefinition,
    ResolvedWorkflow,
    StageDefinition,
)
from txsuite.project.schema import (
    ID_RE,
    SCHEMA_VERSION,
    SUPPORTED_EXECUTION_PROFILES,
    SUPPORTED_MODALITIES,
    WorkflowConfigError,
    parse_artifact_reference,
)


_TOP_LEVEL_KEYS = frozenset({"schema_version", "project", "execution", "workflow"})
_PROJECT_KEYS = frozenset({"id", "modality", "output_root"})
_EXECUTION_KEYS = frozenset({"profile", "resume"})
_WORKFLOW_KEYS = frozenset({"stages"})
_STAGE_KEYS = frozenset(
    {"id", "uses", "depends_on", "inputs", "params", "outputs"}
)


def _reject_unknown_keys(
    value: Mapping[str, Any], allowed: frozenset[str], location: str
) -> None:
    unknown = sorted(str(key) for key in value if key not in allowed)
    if unknown:
        suffix = "s" if len(unknown) != 1 else ""
        raise WorkflowConfigError(
            f"{location} contains unknown key{suffix}: {', '.join(unknown)}"
        )


def _require_table(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise WorkflowConfigError(f"{location} must be a table")
    if any(not isinstance(key, str) for key in value):
        raise WorkflowConfigError(f"{location} keys must be strings")
    return value


def _required(table: Mapping[str, Any], key: str, location: str) -> Any:
    if key not in table:
        raise WorkflowConfigError(f"{location}.{key} is required")
    return table[key]


def _non_empty_string(value: Any, location: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise WorkflowConfigError(f"{location} must be a non-empty string")
    return value


def _id(value: Any, location: str) -> str:
    result = _non_empty_string(value, location)
    if ID_RE.fullmatch(result) is None:
        raise WorkflowConfigError(
            f"{location} must start with a letter and contain only letters, "
            "digits, '_' or '-'"
        )
    return result


def _resolve_path(value: Any, source_dir: Path, location: str) -> Path:
    raw = _non_empty_string(value, location)
    try:
        path = Path(raw)
        if not path.is_absolute():
            path = source_dir / path
        return path.resolve(strict=False)
    except (OSError, RuntimeError, ValueError) as exc:
        raise WorkflowConfigError(f"{location} is not a valid path: {exc}") from exc


def _input_value(value: Any, source_dir: Path, location: str) -> Any:
    if isinstance(value, str):
        reference = parse_artifact_reference(value, location=location)
        return reference or _resolve_path(value, source_dir, location)
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise WorkflowConfigError(f"{location} keys must be strings")
        return {
            key: _input_value(item, source_dir, f"{location}.{key}")
            for key, item in value.items()
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return tuple(
            _input_value(item, source_dir, f"{location}[{index}]")
            for index, item in enumerate(value)
        )
    return _literal_scalar(value, location)


def _literal_scalar(value: Any, location: str) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise WorkflowConfigError(f"{location} must be a finite number")
        return value
    raise WorkflowConfigError(
        f"{location} has unsupported value type {type(value).__name__}; "
        "workflow values must be JSON-serializable"
    )


def _literal_value(value: Any, location: str) -> Any:
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise WorkflowConfigError(f"{location} keys must be strings")
        return {
            key: _literal_value(item, f"{location}.{key}")
            for key, item in value.items()
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return tuple(
            _literal_value(item, f"{location}[{index}]")
            for index, item in enumerate(value)
        )
    return _literal_scalar(value, location)


def _parse_stage(
    raw_stage: Any, index: int, source_dir: Path
) -> StageDefinition:
    location = f"workflow.stages[{index}]"
    stage = _require_table(raw_stage, location)
    _reject_unknown_keys(stage, _STAGE_KEYS, location)

    stage_id = _id(_required(stage, "id", location), f"{location}.id")
    uses = _non_empty_string(_required(stage, "uses", location), f"{location}.uses")

    raw_dependencies = stage.get("depends_on", [])
    if not isinstance(raw_dependencies, list):
        raise WorkflowConfigError(
            f"{location}.depends_on must be an array of stage IDs"
        )
    dependencies = tuple(
        _id(value, f"{location}.depends_on[{dependency_index}]")
        for dependency_index, value in enumerate(raw_dependencies)
    )
    if len(set(dependencies)) != len(dependencies):
        raise WorkflowConfigError(f"{location}.depends_on contains duplicate stage IDs")

    raw_inputs = _require_table(stage.get("inputs", {}), f"{location}.inputs")
    inputs = {
        key: _input_value(value, source_dir, f"{location}.inputs.{key}")
        for key, value in raw_inputs.items()
    }

    raw_params = _require_table(stage.get("params", {}), f"{location}.params")
    params = {
        key: _literal_value(value, f"{location}.params.{key}")
        for key, value in raw_params.items()
    }

    raw_outputs = _require_table(stage.get("outputs", {}), f"{location}.outputs")
    outputs: dict[str, Path] = {}
    for artifact_id, value in raw_outputs.items():
        artifact_id = _id(artifact_id, f"{location}.outputs artifact ID")
        outputs[artifact_id] = _resolve_path(
            value, source_dir, f"{location}.outputs.{artifact_id}"
        )

    return StageDefinition(
        id=stage_id,
        uses=uses,
        depends_on=dependencies,
        inputs=inputs,
        params=params,
        outputs=outputs,
    )


def _walk_references(value: Any):
    if isinstance(value, ArtifactReference):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _walk_references(item)
    elif isinstance(value, tuple):
        for item in value:
            yield from _walk_references(item)


def _validate_graph(stages: tuple[StageDefinition, ...]) -> None:
    by_id: dict[str, StageDefinition] = {}
    for stage in stages:
        if stage.id in by_id:
            raise WorkflowConfigError(f"duplicate workflow stage ID: {stage.id!r}")
        by_id[stage.id] = stage

    for stage in stages:
        for dependency in stage.depends_on:
            if dependency == stage.id:
                raise WorkflowConfigError(f"stage {stage.id!r} cannot depend on itself")
            if dependency not in by_id:
                raise WorkflowConfigError(
                    f"stage {stage.id!r} depends on unknown stage {dependency!r}"
                )
        for reference in _walk_references(stage.inputs):
            if reference.stage_id not in by_id:
                raise WorkflowConfigError(
                    f"stage {stage.id!r} references unknown stage "
                    f"{reference.stage_id!r}"
                )

    # Artifact references are dependency edges even when depends_on omits them.
    # Artifact names are validated later against the selected StageSpec; the
    # outputs table here contains optional path overrides, not declarations.
    state: dict[str, int] = {}

    def visit(stage_id: str, trail: tuple[str, ...]) -> None:
        marker = state.get(stage_id, 0)
        if marker == 2:
            return
        if marker == 1:
            cycle = " -> ".join((*trail, stage_id))
            raise WorkflowConfigError(f"workflow dependency cycle: {cycle}")
        state[stage_id] = 1
        for dependency in by_id[stage_id].dependencies:
            visit(dependency, (*trail, stage_id))
        state[stage_id] = 2

    for stage in stages:
        visit(stage.id, ())


def parse_project_config(
    document: Mapping[str, Any], *, source_path: str | Path
) -> ResolvedWorkflow:
    """Validate a decoded TOML mapping and resolve it into immutable models."""

    root = _require_table(document, "configuration")
    _reject_unknown_keys(root, _TOP_LEVEL_KEYS, "configuration")

    version = _required(root, "schema_version", "configuration")
    if type(version) is not int or version != SCHEMA_VERSION:
        raise WorkflowConfigError(
            f"schema_version must be integer {SCHEMA_VERSION}; got {version!r}"
        )

    resolved_source = Path(source_path).expanduser().resolve(strict=False)
    source_dir = resolved_source.parent

    project_table = _require_table(
        _required(root, "project", "configuration"), "project"
    )
    _reject_unknown_keys(project_table, _PROJECT_KEYS, "project")
    project_id = _id(_required(project_table, "id", "project"), "project.id")
    modality = _non_empty_string(
        _required(project_table, "modality", "project"), "project.modality"
    )
    if modality not in SUPPORTED_MODALITIES:
        choices = ", ".join(sorted(SUPPORTED_MODALITIES))
        raise WorkflowConfigError(
            f"project.modality must be one of {choices}; got {modality!r}"
        )
    project = ProjectDefinition(
        id=project_id,
        modality=modality,
        output_root=_resolve_path(
            _required(project_table, "output_root", "project"),
            source_dir,
            "project.output_root",
        ),
    )

    execution_table = _require_table(
        _required(root, "execution", "configuration"), "execution"
    )
    _reject_unknown_keys(execution_table, _EXECUTION_KEYS, "execution")
    profile = _non_empty_string(
        _required(execution_table, "profile", "execution"), "execution.profile"
    )
    if profile not in SUPPORTED_EXECUTION_PROFILES:
        choices = ", ".join(sorted(SUPPORTED_EXECUTION_PROFILES))
        raise WorkflowConfigError(
            f"execution.profile must be one of {choices}; got {profile!r}"
        )
    resume = _required(execution_table, "resume", "execution")
    if not isinstance(resume, bool):
        raise WorkflowConfigError("execution.resume must be a boolean")
    execution = ExecutionSettings(profile=profile, resume=resume)

    workflow_table = _require_table(
        _required(root, "workflow", "configuration"), "workflow"
    )
    _reject_unknown_keys(workflow_table, _WORKFLOW_KEYS, "workflow")
    raw_stages = _required(workflow_table, "stages", "workflow")
    if not isinstance(raw_stages, list) or not raw_stages:
        raise WorkflowConfigError("workflow.stages must be a non-empty array of tables")
    stages = tuple(
        _parse_stage(raw_stage, index, source_dir)
        for index, raw_stage in enumerate(raw_stages)
    )
    _validate_graph(stages)

    return ResolvedWorkflow(
        schema_version=SCHEMA_VERSION,
        project=project,
        execution=execution,
        stages=stages,
        source_path=resolved_source,
        source_dir=source_dir,
    )


def load_project_config(path: str | Path) -> ResolvedWorkflow:
    """Read, validate, and resolve one canonical v1 project TOML file."""

    source = Path(path).expanduser().resolve(strict=False)
    try:
        with source.open("rb") as handle:
            document = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise WorkflowConfigError(
            f"cannot read project configuration {source}: {exc}"
        ) from exc
    return parse_project_config(document, source_path=source)


# Workflow-oriented aliases for API discoverability.
load_workflow = load_project_config
resolve_project_config = parse_project_config


__all__ = [
    "load_project_config",
    "load_workflow",
    "parse_project_config",
    "resolve_project_config",
]
