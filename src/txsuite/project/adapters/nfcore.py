"""Pinned-release artifact adapters for nf-core result directories.

nf-core output layouts are versioned external interfaces.  Resolvers therefore select
an exact adapter/release pair, never guess a release, and fail when a contract matches
zero or multiple paths.  Every successful resolution records the adapter and release
used as provenance evidence.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Literal

from txsuite.runtime import TxSuiteError


ArtifactKind = Literal["file", "directory", "any"]


@dataclass(frozen=True)
class ArtifactRule:
    """Release-specific relative patterns for one named artifact type."""

    patterns: tuple[str, ...]
    kind: ArtifactKind = "file"

    def __post_init__(self) -> None:
        object.__setattr__(self, "patterns", tuple(self.patterns))
        if not self.patterns:
            raise ValueError("An nf-core artifact rule needs at least one pattern")
        if self.kind not in {"file", "directory", "any"}:
            raise ValueError(f"Unknown nf-core artifact kind: {self.kind!r}")
        for pattern in self.patterns:
            path = PurePosixPath(pattern)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("nf-core artifact patterns must stay below the results root")


@dataclass(frozen=True)
class NfCoreAdapter:
    """Artifact layout adapter bound to one exact nf-core pipeline release."""

    name: str
    pipeline: str
    release: str
    artifacts: Mapping[str, ArtifactRule]

    def __post_init__(self) -> None:
        if not self.name or not self.pipeline or not self.release:
            raise ValueError("nf-core adapter name, pipeline, and release cannot be empty")
        normalized: dict[str, ArtifactRule] = {}
        for artifact_type, rule in self.artifacts.items():
            if not artifact_type or any(char in artifact_type for char in "*?[]/\\"):
                raise ValueError(f"Invalid named artifact type: {artifact_type!r}")
            if not isinstance(rule, ArtifactRule):
                rule = ArtifactRule(tuple(rule))
            normalized[artifact_type] = rule
        object.__setattr__(self, "artifacts", MappingProxyType(normalized))


@dataclass(frozen=True)
class ResolvedNfCoreArtifact:
    """A unique artifact path plus the layout evidence that selected it."""

    artifact_type: str
    path: Path
    adapter: str
    pipeline: str
    release: str
    source: Literal["adapter", "override"]
    pattern: str | None

    @property
    def evidence(self) -> dict[str, str]:
        evidence = {
            "adapter": self.adapter,
            "pipeline": self.pipeline,
            "release": self.release,
            "source": self.source,
        }
        if self.pattern is not None:
            evidence["pattern"] = self.pattern
        return evidence

    def __fspath__(self) -> str:
        return str(self.path)


RNASEQ_3_26_0 = NfCoreAdapter(
    name="nfcore-rnaseq-3.26.0",
    pipeline="nf-core/rnaseq",
    release="3.26.0",
    artifacts={
        "bulk.rnaseq-results": ArtifactRule((".",), "directory"),
        "bulk.gene-counts": ArtifactRule(
            (
                "salmon/salmon.merged.gene_counts.tsv",
                "star_salmon/salmon.merged.gene_counts.tsv",
                "**/salmon.merged.gene_counts.tsv",
                "**/rsem.merged.gene_counts.tsv",
            )
        ),
        "qc.multiqc-report": ArtifactRule(
            ("multiqc/multiqc_report.html", "**/multiqc_report.html")
        ),
    },
)

SCRNASEQ_4_2_0 = NfCoreAdapter(
    name="nfcore-scrnaseq-4.2.0",
    pipeline="nf-core/scrnaseq",
    release="4.2.0",
    artifacts={
        "single-cell.scrnaseq-results": ArtifactRule((".",), "directory"),
        "single-cell.matrix": ArtifactRule(
            ("*/mtx_conversions/combined_matrix.h5ad",),
            "file",
        ),
        "qc.multiqc-report": ArtifactRule(
            ("multiqc/multiqc_report.html", "**/multiqc_report.html")
        ),
    },
)

NFCORE_ADAPTERS: Mapping[str, NfCoreAdapter] = MappingProxyType(
    {adapter.name: adapter for adapter in (RNASEQ_3_26_0, SCRNASEQ_4_2_0)}
)


def list_nfcore_adapters(pipeline: str | None = None) -> tuple[NfCoreAdapter, ...]:
    """List known pinned adapters in deterministic registration order."""

    adapters = tuple(NFCORE_ADAPTERS.values())
    if pipeline is None:
        return adapters
    return tuple(adapter for adapter in adapters if adapter.pipeline == pipeline)


def get_nfcore_adapter(
    *,
    adapter: str | NfCoreAdapter | None = None,
    pipeline: str | None = None,
    release: str | None = None,
) -> NfCoreAdapter:
    """Select one exact adapter and verify any supplied pipeline/release evidence."""

    if isinstance(adapter, NfCoreAdapter):
        selected = adapter
    elif isinstance(adapter, str):
        try:
            selected = NFCORE_ADAPTERS[adapter]
        except KeyError as exc:
            raise TxSuiteError(f"Unknown nf-core artifact adapter {adapter!r}") from exc
    else:
        if pipeline is None or release is None:
            raise TxSuiteError(
                "nf-core artifact resolution requires an adapter or exact pipeline and release"
            )
        matches = tuple(
            candidate
            for candidate in NFCORE_ADAPTERS.values()
            if candidate.pipeline == pipeline and candidate.release == release
        )
        if not matches:
            raise TxSuiteError(
                f"No nf-core artifact adapter for pipeline {pipeline!r} release {release!r}"
            )
        if len(matches) > 1:
            names = ", ".join(candidate.name for candidate in matches)
            raise TxSuiteError(
                f"Ambiguous nf-core adapters for {pipeline!r} {release!r}: {names}"
            )
        selected = matches[0]
    if pipeline is not None and selected.pipeline != pipeline:
        raise TxSuiteError(
            f"Adapter {selected.name!r} is for {selected.pipeline!r}, not {pipeline!r}"
        )
    if release is not None and selected.release != release:
        raise TxSuiteError(
            f"Adapter {selected.name!r} is for release {selected.release!r}, not {release!r}"
        )
    return selected


def _matches_kind(path: Path, kind: ArtifactKind) -> bool:
    if kind == "file":
        return path.is_file()
    if kind == "directory":
        return path.is_dir()
    return path.exists()


def _confined_path(results_dir: Path, path: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(results_dir):
        raise TxSuiteError(f"{label} escapes nf-core results directory {results_dir}")
    return resolved


def _override_path(results_dir: Path, override: str | Path) -> Path:
    path = Path(override).expanduser()
    if not path.is_absolute():
        if ".." in path.parts:
            raise TxSuiteError("Relative nf-core artifact overrides cannot contain '..'")
        path = results_dir / path
    return _confined_path(results_dir, path, "Explicit artifact override")


def _safe_display_path(path: Path, results_dir: Path) -> str:
    try:
        return str(path.relative_to(results_dir))
    except ValueError:
        return str(path)


def resolve_nfcore_artifact(
    results_dir: str | Path,
    artifact_type: str,
    *,
    adapter: str | NfCoreAdapter | None = None,
    pipeline: str | None = None,
    release: str | None = None,
    override: str | Path | None = None,
) -> ResolvedNfCoreArtifact:
    """Resolve exactly one named artifact for a pinned nf-core release."""

    selected = get_nfcore_adapter(adapter=adapter, pipeline=pipeline, release=release)
    try:
        rule = selected.artifacts[artifact_type]
    except KeyError as exc:
        known = ", ".join(selected.artifacts)
        raise TxSuiteError(
            f"Adapter {selected.name!r} does not declare artifact {artifact_type!r}; "
            f"known artifacts: {known}"
        ) from exc
    root = Path(results_dir).expanduser().resolve()
    if not root.is_dir():
        raise TxSuiteError(f"nf-core results directory does not exist: {root}")
    if override is not None:
        path = _override_path(root, override)
        if not _matches_kind(path, rule.kind):
            raise TxSuiteError(
                f"Explicit override for {artifact_type!r} is not an existing "
                f"{rule.kind}: {path}"
            )
        return ResolvedNfCoreArtifact(
            artifact_type=artifact_type,
            path=path,
            adapter=selected.name,
            pipeline=selected.pipeline,
            release=selected.release,
            source="override",
            pattern=None,
        )

    candidates: dict[Path, str] = {}
    for pattern in rule.patterns:
        paths = (root,) if pattern == "." else root.glob(pattern)
        for path in paths:
            if _matches_kind(path, rule.kind):
                confined = _confined_path(root, path, "Resolved artifact")
                candidates.setdefault(confined, pattern)
    ordered = sorted(candidates, key=lambda path: path.as_posix())
    if not ordered:
        raise TxSuiteError(
            f"No {artifact_type!r} artifact found below {root} using adapter "
            f"{selected.name!r} release {selected.release!r}"
        )
    if len(ordered) > 1:
        rendered = ", ".join(_safe_display_path(path, root) for path in ordered)
        raise TxSuiteError(
            f"Ambiguous {artifact_type!r} artifacts for adapter {selected.name!r} "
            f"release {selected.release!r}: {rendered}; provide an explicit override"
        )
    path = ordered[0]
    return ResolvedNfCoreArtifact(
        artifact_type=artifact_type,
        path=path,
        adapter=selected.name,
        pipeline=selected.pipeline,
        release=selected.release,
        source="adapter",
        pattern=candidates[path],
    )


def resolve_nfcore_artifacts(
    results_dir: str | Path,
    *,
    adapter: str | NfCoreAdapter | None = None,
    pipeline: str | None = None,
    release: str | None = None,
    artifact_types: Iterable[str] | None = None,
    overrides: Mapping[str, str | Path] | None = None,
) -> dict[str, ResolvedNfCoreArtifact]:
    """Resolve a requested set (or all) of an adapter's artifact contracts."""

    selected = get_nfcore_adapter(adapter=adapter, pipeline=pipeline, release=release)
    requested = tuple(selected.artifacts if artifact_types is None else artifact_types)
    if len(requested) != len(set(requested)):
        raise TxSuiteError("nf-core artifact_types contains duplicates")
    overrides = {} if overrides is None else dict(overrides)
    unknown_overrides = sorted(set(overrides) - set(requested))
    if unknown_overrides:
        raise TxSuiteError(
            "Overrides supplied for unrequested artifact type(s): "
            + ", ".join(unknown_overrides)
        )
    return {
        artifact_type: resolve_nfcore_artifact(
            results_dir,
            artifact_type,
            adapter=selected,
            override=overrides.get(artifact_type),
        )
        for artifact_type in requested
    }


# Short aliases for integration code that already lives in the nf-core adapter module.
resolve_artifact = resolve_nfcore_artifact
resolve_artifacts = resolve_nfcore_artifacts
NfCoreArtifactResolution = ResolvedNfCoreArtifact


__all__ = [
    "ArtifactKind",
    "ArtifactRule",
    "NFCORE_ADAPTERS",
    "NfCoreAdapter",
    "NfCoreArtifactResolution",
    "RNASEQ_3_26_0",
    "ResolvedNfCoreArtifact",
    "SCRNASEQ_4_2_0",
    "get_nfcore_adapter",
    "list_nfcore_adapters",
    "resolve_artifact",
    "resolve_artifacts",
    "resolve_nfcore_artifact",
    "resolve_nfcore_artifacts",
]
