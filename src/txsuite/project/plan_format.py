"""Stable JSON and concise human renderers for project command plans."""

from __future__ import annotations

import shlex
from collections.abc import Mapping
from typing import Any

from .planner import PlannedWorkflow


def format_plan_json(plan: PlannedWorkflow, *, indent: int | None = 2) -> str:
    """Render a plan as deterministic canonical JSON."""

    return plan.to_json(indent=indent)


def _artifact_lines(label: str, artifacts: Mapping[str, Any]) -> str:
    values: list[str] = []
    for name, artifact in artifacts.items():
        path = artifact.path
        if path is None:
            path = getattr(artifact, "value", getattr(artifact, "reference", "?"))
        values.append(f"{name}<{artifact.artifact_type}>={path}")
    return f"  {label}: " + ("; ".join(values) if values else "-")


def format_plan_human(plan: PlannedWorkflow) -> str:
    """Render the essential DAG, requirements, artifacts, and pins."""

    lines = [
        f"Project {plan.project_id} ({plan.modality}) — {len(plan.stages)} stage(s)",
        f"Output: {plan.output_root}",
    ]
    for index, stage in enumerate(plan.stages, start=1):
        dependencies = ", ".join(stage.depends_on) or "-"
        executables = ", ".join(stage.executables) or "-"
        images = ", ".join(
            f"{key}={value}" for key, value in stage.images.items()
        ) or "-"
        pins = ", ".join(f"{key}={value}" for key, value in stage.pins.items()) or "-"
        postflight = ", ".join(
            f"{key}={value}"
            for key, value in stage.postflight.items()
            if key in {"adapter", "pipeline", "release", "results_root"}
        ) or "-"
        lines.extend(
            [
                "",
                f"{index}. {stage.id} [{stage.uses}; {stage.maturity}; {stage.command_state}]",
                f"  depends: {dependencies}",
                f"  executables: {executables}",
                f"  images: {images}",
                _artifact_lines("inputs", stage.inputs),
                _artifact_lines("outputs", stage.outputs),
                f"  pins: {pins}",
                f"  postflight: {postflight}",
                f"  command: {shlex.join(stage.command)}",
            ]
        )
    return "\n".join(lines)


def render_plan(plan: PlannedWorkflow, *, format: str = "human") -> str:
    """Render ``plan`` as ``human`` or ``json``."""

    if format == "human":
        return format_plan_human(plan)
    if format == "json":
        return format_plan_json(plan)
    raise ValueError("format must be 'human' or 'json'")


# Short aliases for callers that have already selected a plan renderer.
format_json = format_plan_json
format_human = format_plan_human


__all__ = [
    "format_human",
    "format_json",
    "format_plan_human",
    "format_plan_json",
    "render_plan",
]
