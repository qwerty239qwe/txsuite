from __future__ import annotations

import copy
import os
import tomllib
from pathlib import Path
from typing import Any


DEFAULT_CONFIG: dict[str, Any] = {
    "execution": {"profile": "docker"},
    "images": {
        "bulk_r": "txsuite/bulk-r:0.1.0",
        "single_cell_python": "txsuite/single-cell-python:0.1.0",
        "spatial_python": "txsuite/spatial-python:0.1.0",
    },
    "pipelines": {
        "bulk": {"name": "nf-core/rnaseq", "release": "3.26.0"},
        "single_cell": {"name": "nf-core/scrnaseq", "release": "4.2.0"},
        "spatial": {"name": "spaceranger", "release": "4.1.0"},
    },
}

DEFAULT_TOML = """[execution]
profile = "docker"

[images]
bulk_r = "txsuite/bulk-r:0.1.0"
single_cell_python = "txsuite/single-cell-python:0.1.0"
spatial_python = "txsuite/spatial-python:0.1.0"

[pipelines.bulk]
name = "nf-core/rnaseq"
release = "3.26.0"

[pipelines.single_cell]
name = "nf-core/scrnaseq"
release = "4.2.0"

[pipelines.spatial]
name = "spaceranger"
release = "4.1.0"
"""


class ConfigError(ValueError):
    pass


def user_config_path() -> Path:
    if appdata := os.environ.get("APPDATA"):
        return Path(appdata) / "txsuite" / "config.toml"
    base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return base / "txsuite" / "config.toml"


def load_config(
    project_path: Path | None = None,
    *,
    user_path: Path | None = None,
) -> dict[str, Any]:
    config = copy.deepcopy(DEFAULT_CONFIG)
    for path in (user_path or user_config_path(), project_path or Path("txsuite.toml")):
        if path.exists():
            try:
                parsed = tomllib.loads(path.read_text(encoding="utf-8"))
            except (OSError, tomllib.TOMLDecodeError) as exc:
                raise ConfigError(f"Cannot read {path}: {exc}") from exc
            _merge(config, parsed)
    _validate(config)
    return config


def _merge(target: dict[str, Any], source: dict[str, Any]) -> None:
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _merge(target[key], value)
        else:
            target[key] = value


def _validate(config: dict[str, Any]) -> None:
    profile = config.get("execution", {}).get("profile")
    if profile not in {"docker", "apptainer"}:
        raise ConfigError("execution.profile must be 'docker' or 'apptainer'")
    pipelines = config.get("pipelines")
    if not isinstance(pipelines, dict):
        raise ConfigError("pipelines must be a table")
    for modality in ("bulk", "single_cell", "spatial"):
        pipeline = pipelines.get(modality)
        if not isinstance(pipeline, dict) or not all(
            isinstance(pipeline.get(field), str) and pipeline[field]
            for field in ("name", "release")
        ):
            raise ConfigError(
                f"pipelines.{modality} requires non-empty name and release"
            )
    for image in ("bulk_r", "single_cell_python", "spatial_python"):
        if (
            not isinstance(config.get("images", {}).get(image), str)
            or not config["images"][image]
        ):
            raise ConfigError(f"images.{image} must be a non-empty image tag")
