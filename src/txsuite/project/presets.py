"""Discover and safely scaffold packaged TxSuite project presets."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from importlib.resources.abc import Traversable
from pathlib import Path, PurePosixPath

from txsuite.runtime import TxSuiteError


_RESOURCE_PACKAGE = "txsuite.resources.project"


@dataclass(frozen=True, slots=True)
class ProjectPreset:
    """Metadata for one immutable project scaffold bundled with TxSuite."""

    name: str
    description: str
    modality: str
    files: tuple[str, ...]

    @property
    def id(self) -> str:
        """Return the stable preset identifier."""

        return self.name


PROJECT_PRESETS: tuple[ProjectPreset, ...] = (
    ProjectPreset(
        name="bulk-rnaseq",
        description="Bulk RNA-seq, differential expression, and enrichment",
        modality="bulk",
        files=(
            "README.md",
            "genesets.gmt",
            "metadata.tsv",
            "rnaseq-params.json",
            "samplesheet.csv",
            "workflow.toml",
        ),
    ),
    ProjectPreset(
        name="bulk-salmon",
        description="Native salmon quantification followed by differential expression",
        modality="bulk",
        files=(
            "README.md",
            "metadata.tsv",
            "samplesheet.csv",
            "workflow.toml",
        ),
    ),
    ProjectPreset(
        name="scrnaseq",
        description="nf-core/scRNA-seq quantification followed by Scanpy analysis",
        modality="single-cell",
        files=(
            "README.md",
            "samplesheet.csv",
            "scrnaseq-params.json",
            "workflow.toml",
        ),
    ),
    ProjectPreset(
        name="scrnaseq-alevin",
        description="Native simpleaf/alevin-fry quantification followed by Scanpy analysis",
        modality="single-cell",
        files=(
            "README.md",
            "samplesheet.csv",
            "workflow.toml",
        ),
    ),
    ProjectPreset(
        name="scrnaseq-pseudobulk",
        description="Pseudobulk differential expression from an annotated AnnData file",
        modality="single-cell",
        files=(
            "README.md",
            "cell-metadata.tsv",
            "workflow.toml",
        ),
    ),
)

_PRESETS_BY_NAME = {preset.name: preset for preset in PROJECT_PRESETS}


def list_project_presets() -> tuple[ProjectPreset, ...]:
    """Return available presets in stable presentation order."""

    return PROJECT_PRESETS


def get_project_preset(name: str) -> ProjectPreset:
    """Return an exact preset by name, failing clearly for unknown names."""

    try:
        return _PRESETS_BY_NAME[name]
    except (KeyError, TypeError) as exc:
        known = ", ".join(_PRESETS_BY_NAME)
        raise TxSuiteError(
            f"Unknown project preset {name!r}; known presets: {known}"
        ) from exc


def _resource_root(preset: ProjectPreset) -> Traversable:
    root = resources.files(_RESOURCE_PACKAGE).joinpath("presets", preset.name)
    if not root.is_dir():
        raise TxSuiteError(
            f"Packaged resources for project preset {preset.name!r} are missing"
        )
    return root


def _resource_files(
    root: Traversable, relative: PurePosixPath = PurePosixPath()
) -> tuple[tuple[PurePosixPath, bytes], ...]:
    discovered: list[tuple[PurePosixPath, bytes]] = []
    for child in sorted(root.iterdir(), key=lambda item: item.name):
        child_relative = relative / child.name
        if child.is_dir():
            discovered.extend(_resource_files(child, child_relative))
        elif child.is_file():
            discovered.append((child_relative, child.read_bytes()))
    return tuple(discovered)


def scaffold_project_preset(name: str, target: str | Path) -> Path:
    """Copy a packaged preset into a new directory without overwriting anything.

    The target directory must not already exist, including when it is empty. This
    conservative contract lets a future CLI expose an explicit force policy without
    weakening the safe library default.
    """

    preset = get_project_preset(name)
    source_files = _resource_files(_resource_root(preset))
    packaged_names = tuple(path.as_posix() for path, _ in source_files)
    if packaged_names != preset.files:
        raise TxSuiteError(
            f"Packaged resources for project preset {name!r} do not match its manifest"
        )

    destination = Path(target).expanduser()
    try:
        destination.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise FileExistsError(
            f"Project scaffold target already exists: {destination}"
        ) from exc

    try:
        for relative, payload in source_files:
            output = destination.joinpath(*relative.parts)
            output.parent.mkdir(parents=True, exist_ok=True)
            with output.open("xb") as handle:
                handle.write(payload)
    except Exception:
        # The destination was created exclusively above and every file write uses
        # exclusive mode, so removing this partial tree cannot erase prior content.
        import shutil

        shutil.rmtree(destination)
        raise
    return destination


# Concise aliases keep the module pleasant for interactive callers while the
# project-prefixed names remain unambiguous when re-exported by an application.
list_presets = list_project_presets
get_preset = get_project_preset
scaffold_preset = scaffold_project_preset


__all__ = [
    "PROJECT_PRESETS",
    "ProjectPreset",
    "get_preset",
    "get_project_preset",
    "list_presets",
    "list_project_presets",
    "scaffold_preset",
    "scaffold_project_preset",
]
