"""Constants and small helpers for the dependency-free workflow schema."""

from __future__ import annotations

import json
import re
from importlib import resources
from typing import Any, Mapping

from txsuite.project.model import ArtifactReference


SCHEMA_VERSION = 1
SUPPORTED_MODALITIES = frozenset({"bulk", "single-cell", "spatial"})
SUPPORTED_EXECUTION_PROFILES = frozenset({"docker", "apptainer"})
REDACTED = "***REDACTED***"

_ID_PATTERN = r"[A-Za-z][A-Za-z0-9_-]*"
ID_RE = re.compile(rf"^{_ID_PATTERN}$")
ARTIFACT_REFERENCE_RE = re.compile(
    rf"^\$\{{(?P<stage>{_ID_PATTERN})\.(?P<artifact>{_ID_PATTERN})\}}$"
)


class WorkflowConfigError(ValueError):
    """A project workflow is unreadable or violates schema version 1."""


# Compatibility names make the exception easy to discover from either module.
ProjectConfigError = WorkflowConfigError
ConfigError = WorkflowConfigError


def parse_artifact_reference(
    value: str, *, location: str = "value"
) -> ArtifactReference | None:
    """Parse an exact ``${stage.artifact}`` input reference.

    Plain strings return ``None``. Strings that appear to attempt interpolation
    but do not match the v1 grammar are rejected rather than treated as paths.
    """

    match = ARTIFACT_REFERENCE_RE.fullmatch(value)
    if match:
        return ArtifactReference(match.group("stage"), match.group("artifact"))
    if "${" in value:
        raise WorkflowConfigError(
            f"{location} has malformed artifact reference {value!r}; "
            "expected '${stage.artifact}'"
        )
    return None


_SENSITIVE_NAMES = frozenset(
    {
        "api_key",
        "apikey",
        "access_key",
        "authorization",
        "auth_token",
        "client_secret",
        "credential",
        "credentials",
        "password",
        "passwd",
        "private_key",
        "secret",
        "secret_key",
        "token",
    }
)
_SENSITIVE_SEGMENTS = frozenset(
    {"credential", "credentials", "password", "passwd", "secret", "token"}
)


def _is_sensitive_key(key: object) -> bool:
    normalized = re.sub(r"[^a-z0-9]+", "_", str(key).casefold()).strip("_")
    if normalized in _SENSITIVE_NAMES:
        return True
    return bool(set(normalized.split("_")) & _SENSITIVE_SEGMENTS)


def redact_sensitive(value: Any, *, replacement: str = REDACTED) -> Any:
    """Return JSON-like data with sensitive mapping values recursively redacted."""

    if isinstance(value, Mapping):
        return {
            str(key): replacement
            if _is_sensitive_key(key)
            else redact_sensitive(item, replacement=replacement)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [redact_sensitive(item, replacement=replacement) for item in value]
    return value


def load_workflow_schema() -> dict[str, Any]:
    """Load a fresh copy of the packaged JSON Schema for workflow v1."""

    schema = resources.files("txsuite.resources.project").joinpath(
        "workflow.schema.json"
    )
    return json.loads(schema.read_text(encoding="utf-8"))


__all__ = [
    "ARTIFACT_REFERENCE_RE",
    "ConfigError",
    "ID_RE",
    "ProjectConfigError",
    "REDACTED",
    "SCHEMA_VERSION",
    "SUPPORTED_EXECUTION_PROFILES",
    "SUPPORTED_MODALITIES",
    "WorkflowConfigError",
    "load_workflow_schema",
    "parse_artifact_reference",
    "redact_sensitive",
]
