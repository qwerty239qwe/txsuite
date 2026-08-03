"""Immutable value objects for a resolved TxSuite project workflow."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping


def _freeze(value: Any) -> Any:
    """Recursively make supported configuration values immutable."""

    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze(item) for key, item in value.items()}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _json_value(value: Any) -> Any:
    """Convert resolved values to fresh JSON-compatible containers."""

    if isinstance(value, ArtifactReference):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _artifact_references(value: Any):
    if isinstance(value, ArtifactReference):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _artifact_references(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _artifact_references(item)


@dataclass(frozen=True, slots=True)
class ArtifactReference:
    """A reference from a stage input to another stage's named output."""

    stage_id: str
    artifact_id: str

    @property
    def stage(self) -> str:
        """Return the referenced stage ID."""

        return self.stage_id

    @property
    def artifact(self) -> str:
        """Return the referenced output artifact ID."""

        return self.artifact_id

    def __str__(self) -> str:
        return f"${{{self.stage_id}.{self.artifact_id}}}"

    def to_dict(self) -> dict[str, str]:
        """Return a structured description useful to API consumers."""

        return {"stage": self.stage_id, "artifact": self.artifact_id}


@dataclass(frozen=True, slots=True)
class ProjectDefinition:
    """Resolved project identity and output location."""

    id: str
    modality: str
    output_root: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "output_root", Path(self.output_root))

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "modality": self.modality,
            "output_root": str(self.output_root),
        }


@dataclass(frozen=True, slots=True)
class ExecutionSettings:
    """Execution settings shared by every stage in a project."""

    profile: str
    resume: bool

    def to_dict(self) -> dict[str, Any]:
        return {"profile": self.profile, "resume": self.resume}


@dataclass(frozen=True, slots=True)
class StageDefinition:
    """One immutable stage in the project workflow DAG."""

    id: str
    uses: str
    depends_on: tuple[str, ...] = ()
    inputs: Mapping[str, Any] = field(default_factory=dict)
    params: Mapping[str, Any] = field(default_factory=dict)
    outputs: Mapping[str, Path] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "depends_on", tuple(self.depends_on))
        object.__setattr__(self, "inputs", _freeze(self.inputs))
        object.__setattr__(self, "params", _freeze(self.params))
        object.__setattr__(
            self,
            "outputs",
            MappingProxyType(
                {str(key): Path(value) for key, value in self.outputs.items()}
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "uses": self.uses,
            "depends_on": list(self.depends_on),
            "inputs": _json_value(self.inputs),
            "params": _json_value(self.params),
            "outputs": _json_value(self.outputs),
        }

    @property
    def artifact_references(self) -> tuple[ArtifactReference, ...]:
        """Return input artifact references in deterministic traversal order."""

        return tuple(_artifact_references(self.inputs))

    @property
    def dependencies(self) -> tuple[str, ...]:
        """Return explicit and artifact-implied dependency stage IDs."""

        ordered = (
            *self.depends_on,
            *(ref.stage_id for ref in self.artifact_references),
        )
        return tuple(dict.fromkeys(ordered))


@dataclass(frozen=True, slots=True)
class ResolvedWorkflow:
    """A validated workflow with all filesystem paths made absolute."""

    schema_version: int
    project: ProjectDefinition
    execution: ExecutionSettings
    stages: tuple[StageDefinition, ...]
    source_path: Path
    source_dir: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "stages", tuple(self.stages))
        object.__setattr__(self, "source_path", Path(self.source_path))
        object.__setattr__(self, "source_dir", Path(self.source_dir))

    def stage(self, stage_id: str) -> StageDefinition:
        """Return a stage by ID, raising ``KeyError`` when it is absent."""

        for stage in self.stages:
            if stage.id == stage_id:
                return stage
        raise KeyError(stage_id)

    def to_dict(self, *, redact: bool = False) -> dict[str, Any]:
        """Return the canonical, JSON-serializable resolved representation."""

        value = {
            "schema_version": self.schema_version,
            "source_path": str(self.source_path),
            "source_dir": str(self.source_dir),
            "project": self.project.to_dict(),
            "execution": self.execution.to_dict(),
            "workflow": {"stages": [stage.to_dict() for stage in self.stages]},
        }
        if redact:
            # Local import avoids a model/schema import cycle.
            from txsuite.project.schema import redact_sensitive

            return redact_sensitive(value)
        return value

    def to_json(self, *, indent: int | None = None, redact: bool = False) -> str:
        """Serialize the canonical representation with deterministic key order."""

        separators = (",", ":") if indent is None else None
        return json.dumps(
            self.to_dict(redact=redact),
            sort_keys=True,
            indent=indent,
            separators=separators,
            allow_nan=False,
        )


# Concise compatibility aliases for callers that prefer domain nouns.
Project = ProjectDefinition
Execution = ExecutionSettings
Stage = StageDefinition


__all__ = [
    "ArtifactReference",
    "Execution",
    "ExecutionSettings",
    "Project",
    "ProjectDefinition",
    "ResolvedWorkflow",
    "Stage",
    "StageDefinition",
]
