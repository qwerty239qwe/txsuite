"""Pure, deterministic planning for resolved TxSuite project workflows."""

from __future__ import annotations

import copy
import heapq
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Literal

from txsuite.runtime import TxSuiteError

from .model import ArtifactReference, ResolvedWorkflow, StageDefinition
from .registry import StageSpec, get_stage_spec, validate_stage_parameters


CommandState = Literal["resolved", "deferred"]
InputSource = Literal["literal", "artifact", "mixed"]


class PlanningError(TxSuiteError):
    """Raised when a resolved workflow cannot form a runnable command plan."""


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType(
            {str(key): _freeze(item) for key, item in sorted(value.items())}
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


def _json_value(value: Any) -> Any:
    if isinstance(value, ArtifactReference):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _references(value: Any):
    if isinstance(value, ArtifactReference):
        yield value
    elif isinstance(value, Mapping):
        for item in value.values():
            yield from _references(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _references(item)


def _contains_literal(value: Any) -> bool:
    if isinstance(value, ArtifactReference):
        return False
    if isinstance(value, Mapping):
        return any(_contains_literal(item) for item in value.values())
    if isinstance(value, (list, tuple)):
        return any(_contains_literal(item) for item in value)
    return True


@dataclass(frozen=True, slots=True)
class PlannedInput:
    """One typed input after all currently knowable references are resolved."""

    artifact_type: str
    value: Any
    source: InputSource = "literal"
    references: tuple[ArtifactReference, ...] = ()
    deferred: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", _freeze(self.value))
        object.__setattr__(self, "references", tuple(self.references))

    @property
    def path(self) -> Path | None:
        return self.value if isinstance(self.value, Path) else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "type": self.artifact_type,
            "source": self.source,
            "value": _json_value(self.value),
            "references": [reference.to_dict() for reference in self.references],
            "deferred": self.deferred,
        }


@dataclass(frozen=True, slots=True)
class PlannedOutput:
    """One named output, concrete when its location is known at plan time."""

    stage_id: str
    name: str
    artifact_type: str
    kind: str = "file"
    non_empty: bool = False
    path: Path | None = None
    explicit: bool = False

    def __post_init__(self) -> None:
        if self.kind not in {"file", "directory"}:
            raise ValueError("PlannedOutput kind must be 'file' or 'directory'")
        if not isinstance(self.non_empty, bool):
            raise TypeError("PlannedOutput non_empty must be a boolean")
        if self.path is not None:
            object.__setattr__(self, "path", Path(self.path))

    @property
    def reference(self) -> ArtifactReference:
        return ArtifactReference(self.stage_id, self.name)

    @property
    def value(self) -> Path | ArtifactReference:
        return self.path if self.path is not None else self.reference

    def to_dict(self) -> dict[str, Any]:
        return {
            "type": self.artifact_type,
            "kind": self.kind,
            "non_empty": self.non_empty,
            "path": None if self.path is None else str(self.path),
            "reference": str(self.reference),
            "explicit": self.explicit,
        }


@dataclass(frozen=True, slots=True)
class PlannedStage:
    """An immutable stage in deterministic execution order."""

    id: str
    uses: str
    modality: str
    maturity: str
    depends_on: tuple[str, ...]
    inputs: Mapping[str, PlannedInput]
    outputs: Mapping[str, PlannedOutput]
    params: Mapping[str, Any]
    outdir: Path
    command: tuple[str, ...]
    command_state: CommandState
    executables: tuple[str, ...] = ()
    images: Mapping[str, str] = field(default_factory=dict)
    pins: Mapping[str, str] = field(default_factory=dict)
    postflight: Mapping[str, Any] = field(default_factory=dict)
    supports_resume: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "depends_on", tuple(self.depends_on))
        object.__setattr__(self, "inputs", _freeze(self.inputs))
        object.__setattr__(self, "outputs", _freeze(self.outputs))
        object.__setattr__(self, "params", _freeze(self.params))
        object.__setattr__(self, "outdir", Path(self.outdir))
        object.__setattr__(self, "command", tuple(self.command))
        object.__setattr__(self, "executables", tuple(self.executables))
        object.__setattr__(self, "images", _freeze(self.images))
        object.__setattr__(self, "pins", _freeze(self.pins))
        object.__setattr__(self, "postflight", _freeze(self.postflight))

    @property
    def resolved_inputs(self) -> Mapping[str, Any]:
        return MappingProxyType(
            {name: planned.value for name, planned in self.inputs.items()}
        )

    @property
    def required_outputs(self) -> tuple[Mapping[str, Any], ...]:
        return tuple(
            MappingProxyType(
                {
                    "path": output.path,
                    "kind": output.kind,
                    "non_empty": output.non_empty,
                }
            )
            for output in self.outputs.values()
            if output.path is not None
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "uses": self.uses,
            "modality": self.modality,
            "maturity": self.maturity,
            "depends_on": list(self.depends_on),
            "inputs": {
                name: _json_value(planned.value)
                for name, planned in self.inputs.items()
            },
            "input_artifacts": {
                name: planned.to_dict() for name, planned in self.inputs.items()
            },
            "params": _json_value(self.params),
            "outdir": str(self.outdir),
            "outputs": {
                name: _json_value(planned.value)
                for name, planned in self.outputs.items()
            },
            "output_artifacts": {
                name: planned.to_dict() for name, planned in self.outputs.items()
            },
            "required_outputs": _json_value(self.required_outputs),
            "command": list(self.command),
            "command_state": self.command_state,
            "requirements": {
                "executables": list(self.executables),
                "images": dict(self.images),
            },
            "pins": dict(self.pins),
            "postflight": _json_value(self.postflight),
            "supports_resume": self.supports_resume,
        }


@dataclass(frozen=True, slots=True)
class PlannedWorkflow:
    """Canonical, immutable command plan for one resolved workflow."""

    schema_version: int
    project_id: str
    modality: str
    output_root: Path
    profile: str
    resume: bool
    source_path: Path
    stages: tuple[PlannedStage, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "output_root", Path(self.output_root))
        object.__setattr__(self, "source_path", Path(self.source_path))
        object.__setattr__(self, "stages", tuple(self.stages))

    def stage(self, stage_id: str) -> PlannedStage:
        for stage in self.stages:
            if stage.id == stage_id:
                return stage
        raise KeyError(stage_id)

    def to_dict(self) -> dict[str, Any]:
        executables = sorted(
            {item for stage in self.stages for item in stage.executables}
        )
        images = sorted(
            {value for stage in self.stages for value in stage.images.values()}
        )
        return {
            "schema_version": self.schema_version,
            "project": {
                "id": self.project_id,
                "modality": self.modality,
                "output_root": str(self.output_root),
            },
            "execution": {"profile": self.profile, "resume": self.resume},
            "source_path": str(self.source_path),
            "requirements": {"executables": executables, "images": images},
            "stages": [stage.to_dict() for stage in self.stages],
        }

    def to_json(self, *, indent: int | None = None) -> str:
        separators = (",", ":") if indent is None else None
        return json.dumps(
            self.to_dict(),
            sort_keys=True,
            indent=indent,
            separators=separators,
            allow_nan=False,
        )

    def to_human(self) -> str:
        from .plan_format import format_plan_human

        return format_plan_human(self)


def _topological_order(workflow: ResolvedWorkflow) -> tuple[StageDefinition, ...]:
    by_id: dict[str, StageDefinition] = {}
    for stage in workflow.stages:
        if stage.id in by_id:
            raise PlanningError(f"Duplicate workflow stage ID {stage.id!r}")
        by_id[stage.id] = stage

    dependencies: dict[str, tuple[str, ...]] = {}
    dependents: dict[str, list[str]] = {stage_id: [] for stage_id in by_id}
    for stage in workflow.stages:
        deps = tuple(dict.fromkeys(stage.dependencies))
        for dependency in deps:
            if dependency not in by_id:
                raise PlanningError(
                    f"Stage {stage.id!r} depends on unknown stage {dependency!r}"
                )
            dependents[dependency].append(stage.id)
        dependencies[stage.id] = deps

    ready = [stage_id for stage_id, deps in dependencies.items() if not deps]
    heapq.heapify(ready)
    ordered: list[StageDefinition] = []
    remaining = {stage_id: len(deps) for stage_id, deps in dependencies.items()}
    while ready:
        stage_id = heapq.heappop(ready)
        ordered.append(by_id[stage_id])
        for consumer in sorted(dependents[stage_id]):
            remaining[consumer] -= 1
            if remaining[consumer] == 0:
                heapq.heappush(ready, consumer)
    if len(ordered) != len(by_id):
        cyclic = ", ".join(sorted(stage_id for stage_id, count in remaining.items() if count))
        raise PlanningError(f"Workflow dependency cycle involves: {cyclic}")
    return tuple(ordered)


def _absolute(path: Path, source_dir: Path) -> Path:
    return path if path.is_absolute() else (source_dir / path).resolve(strict=False)


def _resolve_parameter_paths(value: Any, source_dir: Path) -> Any:
    """Anchor validator-produced paths to the workflow, preserving strings."""

    if isinstance(value, Path):
        return _absolute(value, source_dir)
    if isinstance(value, Mapping):
        return {
            str(key): _resolve_parameter_paths(item, source_dir)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return tuple(_resolve_parameter_paths(item, source_dir) for item in value)
    return value


def _default_output_path(
    uses: str,
    name: str,
    params: Mapping[str, Any],
    outdir: Path,
    inputs: Mapping[str, Any],
) -> Path | None:
    if name == "results":
        return outdir
    known = {
        # The native quantification DAGs publish fixed relative paths, so unlike
        # the nf-core stages their artifacts resolve at plan time.
        ("bulk.salmon", "counts"): outdir / "counts" / "gene_counts.tsv",
        ("bulk.salmon", "tx_counts"): outdir / "counts" / "transcript_counts.tsv",
        ("single-cell.alevin", "matrix"): outdir / "matrix" / "alevin.h5ad",
        ("bulk.align", "counts"): outdir / "counts" / "gene_counts.tsv",
        ("bulk.align", "strandedness"): outdir / "counts" / "strandedness.tsv",
        ("bulk.align", "alignments"): outdir / "alignments",
        ("bulk.align", "logs"): outdir / "logs",
        ("bulk.star-reference", "star_index"): outdir / "star_index",
        ("bulk.star-reference", "manifest"): outdir
        / "reference"
        / "reference-manifest.tsv",
        ("single-cell.scanpy", "h5ad"): outdir / "analysis.h5ad",
        ("single-cell.pseudobulk", "counts"): outdir / "pseudobulk-counts.tsv",
        ("single-cell.pseudobulk", "metadata"): outdir / "pseudobulk-metadata.tsv",
        ("single-cell.pseudobulk-de", "de_results"): outdir
        / "deseq2"
        / "deseq2-results.tsv",
    }
    if uses == "bulk.star-reference" and name in {"fasta_fai", "dict"}:
        # GATK looks for <fasta>.fai and <fasta minus extension>.dict beside the
        # reference, so these names follow the caller's FASTA rather than a
        # canonical one. A FASTA produced by another stage has no name at plan
        # time, which the caller reports rather than silently deferring.
        fasta = inputs.get("fasta")
        if not isinstance(fasta, Path):
            return None
        suffix = f"{fasta.name}.fai" if name == "fasta_fai" else f"{fasta.stem}.dict"
        return outdir / "reference" / suffix
    if uses == "bulk.de" and name == "de_results":
        return outdir / f"{params['method']}-results.tsv"
    if uses == "bulk.de" and name == "contrast_index":
        return outdir / "contrasts.tsv"
    return known.get((uses, name))


def _resolve_outputs(
    stage: StageDefinition,
    spec: StageSpec,
    params: Mapping[str, Any],
    output_root: Path,
    source_dir: Path,
) -> tuple[Path, Mapping[str, PlannedOutput]]:
    unknown = sorted(set(stage.outputs) - set(spec.outputs))
    if unknown:
        raise PlanningError(
            f"Stage {stage.id!r} overrides unknown output(s) for {spec.uses}: "
            + ", ".join(unknown)
        )
    default_outdir = output_root / stage.id
    outdir = _absolute(stage.outputs.get("results", default_outdir), source_dir)
    outputs: dict[str, PlannedOutput] = {}
    for name, artifact_type in sorted(spec.outputs.items()):
        explicit = name in stage.outputs
        path = (
            _absolute(stage.outputs[name], source_dir)
            if explicit
            else _default_output_path(spec.uses, name, params, outdir, stage.inputs)
        )
        outputs[name] = PlannedOutput(
            stage_id=stage.id,
            name=name,
            artifact_type=artifact_type,
            kind=spec.output_policies[name].kind,
            non_empty=spec.output_policies[name].non_empty,
            path=path,
            explicit=explicit,
        )
    return outdir, MappingProxyType(outputs)


def _replace_references(
    value: Any, outputs: Mapping[str, Mapping[str, PlannedOutput]]
) -> tuple[Any, bool]:
    if isinstance(value, ArtifactReference):
        output = outputs[value.stage_id][value.artifact_id]
        return output.value, output.path is None
    if isinstance(value, Mapping):
        resolved: dict[str, Any] = {}
        deferred = False
        for key, item in value.items():
            resolved_item, item_deferred = _replace_references(item, outputs)
            resolved[str(key)] = resolved_item
            deferred = deferred or item_deferred
        return resolved, deferred
    if isinstance(value, (list, tuple)):
        resolved_items = [_replace_references(item, outputs) for item in value]
        return tuple(item for item, _ in resolved_items), any(
            deferred for _, deferred in resolved_items
        )
    return value, False


def _resolve_inputs(
    stage: StageDefinition,
    spec: StageSpec,
    outputs: Mapping[str, Mapping[str, PlannedOutput]],
) -> Mapping[str, PlannedInput]:
    accepted = spec.accepted_inputs
    missing = sorted(set(spec.inputs) - set(stage.inputs))
    unknown = sorted(set(stage.inputs) - set(accepted))
    if missing:
        raise PlanningError(
            f"Stage {stage.id!r} is missing required input(s): {', '.join(missing)}"
        )
    if unknown:
        raise PlanningError(
            f"Stage {stage.id!r} has unknown input(s) for {spec.uses}: "
            + ", ".join(unknown)
        )
    planned: dict[str, PlannedInput] = {}
    # Optional inputs are planned exactly like required ones when supplied, and
    # simply absent otherwise, so a stage that can reuse a prebuilt artifact
    # keeps full reference resolution without forcing every project to provide it.
    for name, artifact_type in sorted(accepted.items()):
        if name not in stage.inputs:
            continue
        raw = stage.inputs[name]
        references = tuple(_references(raw))
        for reference in references:
            producer = outputs.get(reference.stage_id)
            if producer is None:
                raise PlanningError(
                    f"Stage {stage.id!r} references unknown stage {reference.stage_id!r}"
                )
            if reference.artifact_id not in producer:
                raise PlanningError(
                    f"Stage {stage.id!r} input {name!r} references unknown artifact "
                    f"{reference.stage_id}.{reference.artifact_id}"
                )
            actual = producer[reference.artifact_id].artifact_type
            if actual != artifact_type:
                raise PlanningError(
                    f"Artifact type mismatch for stage {stage.id!r} input {name!r}: "
                    f"expected {artifact_type!r}, got {actual!r} from "
                    f"{reference.stage_id}.{reference.artifact_id}"
                )
        value, deferred = _replace_references(raw, outputs)
        source: InputSource
        if not references:
            source = "literal"
        elif _contains_literal(raw):
            source = "mixed"
        else:
            source = "artifact"
        planned[name] = PlannedInput(
            artifact_type=artifact_type,
            value=value,
            source=source,
            references=references,
            deferred=deferred,
        )
    return MappingProxyType(planned)


# Stages that launch a pinned upstream pipeline, and the config key holding its
# pin. Kept in one place so a new launcher cannot pick up pins without also
# getting postflight artifact resolution.
_PIPELINE_KEYS = {
    "bulk.rnaseq": "pipelines.bulk",
    "bulk.rnavar": "pipelines.variants",
    "single-cell.scrnaseq": "pipelines.single_cell",
}


def _config_value(config: Mapping[str, Any], dotted: str) -> Any:
    value: Any = config
    for part in dotted.split("."):
        if not isinstance(value, Mapping) or part not in value:
            raise PlanningError(f"Global config is missing {dotted}")
        value = value[part]
    return value


def _requirements(
    spec: StageSpec, params: Mapping[str, Any], config: Mapping[str, Any]
) -> tuple[Mapping[str, str], Mapping[str, str]]:
    images: dict[str, str] = {}
    for key in spec.required_images:
        value = _config_value(config, key)
        override_name = {
            "images.bulk_r": "bulk_image",
            "images.single_cell_python": "single_cell_image",
        }.get(key)
        if len(spec.required_images) == 1 and params.get("image") is not None:
            value = params["image"]
        elif override_name is not None and params.get(override_name) is not None:
            value = params[override_name]
        if not isinstance(value, str) or not value.strip():
            raise PlanningError(f"Image pin {key} must be a non-empty string")
        images[key] = value

    pins = dict(images)
    pipeline_key = _PIPELINE_KEYS.get(spec.uses)
    if pipeline_key is not None:
        for field_name in ("name", "release"):
            dotted = f"{pipeline_key}.{field_name}"
            value = _config_value(config, dotted)
            if not isinstance(value, str) or not value.strip():
                raise PlanningError(f"Pipeline pin {dotted} must be a non-empty string")
            pins[dotted] = value
    return MappingProxyType(images), MappingProxyType(pins)


def _command_config(
    config: Mapping[str, Any], workflow: ResolvedWorkflow
) -> dict[str, Any]:
    # A stage builder may consume ordinary mutable dictionaries, but planning must
    # neither mutate nor retain aliases into the caller's merged configuration.
    def clone(value: Any) -> Any:
        if isinstance(value, Mapping):
            return {str(key): clone(item) for key, item in value.items()}
        if isinstance(value, list):
            return [clone(item) for item in value]
        if isinstance(value, tuple):
            return tuple(clone(item) for item in value)
        return copy.deepcopy(value)

    resolved = clone(config)
    execution = dict(resolved.get("execution", {}))
    execution["profile"] = workflow.execution.profile
    resolved["execution"] = execution
    return resolved


def _validate_profile(stage: StageDefinition, spec: StageSpec, profile: str) -> None:
    if "nextflow" in spec.required_executables:
        # Nextflow selects the container runtime through its own profile, so any
        # stage that shells out to Nextflow works under both.
        supported = ("apptainer", "docker")
    else:
        # The owned downstream adapters construct `docker run` argv directly.
        supported = ("docker",)
    if profile not in supported:
        choices = ", ".join(supported)
        detail = (
            "; this owned downstream adapter currently emits Docker commands"
            if supported == ("docker",)
            else ""
        )
        raise PlanningError(
            f"Stage {stage.id!r} ({spec.uses}) does not support execution profile "
            f"{profile!r}; supported profile(s): {choices}{detail}"
        )


def _build_command(
    spec: StageSpec,
    inputs: Mapping[str, PlannedInput],
    params: Mapping[str, Any],
    outdir: Path,
    workflow: ResolvedWorkflow,
    config: Mapping[str, Any],
) -> tuple[tuple[str, ...], CommandState]:
    deferred = any(planned.deferred for planned in inputs.values())

    def command_value(value: Any) -> Any:
        if isinstance(value, ArtifactReference):
            return Path(str(value))
        if isinstance(value, Mapping):
            return {str(key): command_value(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return tuple(command_value(item) for item in value)
        return value

    command_inputs = {
        name: command_value(planned.value) for name, planned in inputs.items()
    }
    symbolic_paths: dict[str, str] = {}
    for planned in inputs.values():
        for reference in _references(planned.value):
            token = str(reference)
            symbolic_paths[str(Path(token).resolve(strict=False))] = token
    try:
        command = spec.build_command(
            {
                "global_config": config,
                "resolved_inputs": command_inputs,
                "params": params,
                "outdir": outdir,
                "resume": workflow.execution.resume and spec.supports_resume,
                # Command adapters honor this planner-only seam without performing
                # input existence checks.  It never weakens execution validation.
                "check_inputs": False,
                "planning": True,
            }
        )
    except (KeyError, TypeError, ValueError, OSError, TxSuiteError) as exc:
        raise PlanningError(
            f"Cannot build command for stage {workflow.project.id}/{spec.uses}: {exc}"
        ) from exc
    normalized: list[str] = []
    for argument in command:
        value = str(argument)
        for concrete, token in sorted(symbolic_paths.items(), key=lambda item: -len(item[0])):
            value = value.replace(concrete, token)
        normalized.append(value)
    return tuple(normalized), "deferred" if deferred else "resolved"


def _postflight_metadata(
    spec: StageSpec,
    pins: Mapping[str, str],
    outputs: Mapping[str, PlannedOutput],
    outdir: Path,
) -> Mapping[str, Any]:
    pipeline_key = _PIPELINE_KEYS.get(spec.uses)
    if pipeline_key is None:
        return MappingProxyType({})
    return _freeze(
        {
            "adapter": "nfcore",
            "pipeline": pins[f"{pipeline_key}.name"],
            "release": pins[f"{pipeline_key}.release"],
            "results_root": outdir,
            "artifacts": {
                name: output.artifact_type for name, output in outputs.items()
            },
            "policies": {
                name: {"kind": output.kind, "non_empty": output.non_empty}
                for name, output in outputs.items()
            },
            "overrides": {
                name: output.path
                for name, output in outputs.items()
                if output.explicit
            },
        }
    )


def plan_workflow(
    workflow: ResolvedWorkflow, global_config: Mapping[str, Any]
) -> PlannedWorkflow:
    """Validate and plan ``workflow`` without writes, downloads, or subprocesses."""

    if not isinstance(workflow, ResolvedWorkflow):
        raise TypeError("workflow must be a ResolvedWorkflow")
    if not isinstance(global_config, Mapping):
        raise TypeError("global_config must be a mapping")

    ordered = _topological_order(workflow)
    effective_config = _command_config(global_config, workflow)
    specs: dict[str, StageSpec] = {}
    params_by_stage: dict[str, Mapping[str, Any]] = {}
    outputs_by_stage: dict[str, Mapping[str, PlannedOutput]] = {}
    outdirs: dict[str, Path] = {}
    for stage in ordered:
        try:
            spec = get_stage_spec(stage.uses)
        except TxSuiteError as exc:
            raise PlanningError(f"Stage {stage.id!r}: {exc}") from exc
        if spec.modality != workflow.project.modality:
            raise PlanningError(
                f"Stage {stage.id!r} uses modality {spec.modality!r}, but project "
                f"modality is {workflow.project.modality!r}"
            )
        _validate_profile(stage, spec, workflow.execution.profile)
        try:
            params = validate_stage_parameters(spec, stage.params)
        except TxSuiteError as exc:
            raise PlanningError(f"Stage {stage.id!r}: {exc}") from exc
        params = _resolve_parameter_paths(params, workflow.source_dir)
        outdir, outputs = _resolve_outputs(
            stage,
            spec,
            params,
            workflow.project.output_root,
            workflow.source_dir,
        )
        specs[stage.id] = spec
        params_by_stage[stage.id] = _freeze(params)
        outdirs[stage.id] = outdir
        outputs_by_stage[stage.id] = outputs

    planned_stages: list[PlannedStage] = []
    for stage in ordered:
        spec = specs[stage.id]
        params = params_by_stage[stage.id]
        inputs = _resolve_inputs(stage, spec, outputs_by_stage)
        images, pins = _requirements(spec, params, effective_config)
        command, state = _build_command(
            spec,
            inputs,
            params,
            outdirs[stage.id],
            workflow,
            effective_config,
        )
        planned_stages.append(
            PlannedStage(
                id=stage.id,
                uses=spec.uses,
                modality=spec.modality,
                maturity=spec.maturity,
                depends_on=tuple(sorted(set(stage.dependencies))),
                inputs=inputs,
                outputs=outputs_by_stage[stage.id],
                params=params,
                outdir=outdirs[stage.id],
                command=command,
                command_state=state,
                executables=spec.required_executables,
                images=images,
                pins=pins,
                postflight=_postflight_metadata(
                    spec, pins, outputs_by_stage[stage.id], outdirs[stage.id]
                ),
                supports_resume=spec.supports_resume,
            )
        )
    return PlannedWorkflow(
        schema_version=workflow.schema_version,
        project_id=workflow.project.id,
        modality=workflow.project.modality,
        output_root=workflow.project.output_root,
        profile=workflow.execution.profile,
        resume=workflow.execution.resume,
        source_path=workflow.source_path,
        stages=tuple(planned_stages),
    )


build_plan = plan_workflow
plan_project = plan_workflow


__all__ = [
    "PlannedInput",
    "PlannedOutput",
    "PlannedStage",
    "PlannedWorkflow",
    "PlanningError",
    "build_plan",
    "plan_project",
    "plan_workflow",
]
