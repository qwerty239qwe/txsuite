"""Command-adapter helpers for project stages.

Adapters intentionally consume a small mapping instead of project-model classes.  This
keeps the stage registry usable by the planner, tests, and third-party callers without
coupling it to model construction details.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError


_MISSING = object()


def context_value(context: Mapping[str, Any] | object, name: str, default: Any = _MISSING) -> Any:
    """Read a value from a mapping context or an attribute-based context."""

    if isinstance(context, Mapping):
        if name in context:
            return context[name]
    elif hasattr(context, name):
        return getattr(context, name)
    if default is not _MISSING:
        return default
    raise TxSuiteError(f"Stage command context is missing {name!r}")


def command_context(
    context: Mapping[str, Any] | object,
) -> tuple[
    Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Path, bool, bool
]:
    """Return config, inputs, params, outdir, resume, and input-check policy.

    ``global_config``/``resolved_inputs`` are the canonical planner names.  The shorter
    ``config``/``inputs`` spellings remain accepted for simple direct callers.
    """

    config = context_value(context, "global_config", None)
    if config is None:
        config = context_value(context, "config")
    inputs = context_value(context, "resolved_inputs", None)
    if inputs is None:
        inputs = context_value(context, "inputs")
    params = context_value(context, "params", {})
    outdir = Path(context_value(context, "outdir"))
    resume = context_value(context, "resume", False)
    check_inputs = context_value(context, "check_inputs", True)
    if not isinstance(config, Mapping):
        raise TxSuiteError("Stage command context config must be a mapping")
    if not isinstance(inputs, Mapping):
        raise TxSuiteError("Stage command context inputs must be a mapping")
    if not isinstance(params, Mapping):
        raise TxSuiteError("Stage command context params must be a mapping")
    if not isinstance(resume, bool):
        raise TxSuiteError("Stage command context resume must be a boolean")
    if not isinstance(check_inputs, bool):
        raise TxSuiteError("Stage command context check_inputs must be a boolean")
    return config, inputs, params, outdir, resume, check_inputs


def input_path(inputs: Mapping[str, Any], name: str) -> Path:
    """Resolve a named input, accepting paths and lightweight artifact records."""

    if name not in inputs:
        raise TxSuiteError(f"Stage command context is missing input {name!r}")
    value = inputs[name]
    if isinstance(value, Mapping) and "path" in value:
        value = value["path"]
    elif not isinstance(value, (str, bytes, Path)) and hasattr(value, "path"):
        value = value.path
    try:
        return Path(value)
    except TypeError as exc:
        raise TxSuiteError(f"Stage input {name!r} does not contain a filesystem path") from exc


def optional_input_path(inputs: Mapping[str, Any], name: str) -> Path | None:
    """Resolve an optional input, returning ``None`` when it was not supplied."""

    if name not in inputs:
        return None
    return input_path(inputs, name)


def configured_image(
    config: Mapping[str, Any], params: Mapping[str, Any], image_key: str
) -> str:
    """Choose an explicit stage image or the configured image and reject blanks."""

    image = params.get("image")
    if image is None:
        images = config.get("images")
        if not isinstance(images, Mapping) or image_key not in images:
            raise TxSuiteError(f"Global config is missing images.{image_key}")
        image = images[image_key]
    if not isinstance(image, str) or not image.strip():
        raise TxSuiteError(f"Container image images.{image_key} cannot be empty")
    return image


__all__ = [
    "command_context",
    "configured_image",
    "context_value",
    "input_path",
    "optional_input_path",
]
