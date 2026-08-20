"""Materialized gene-set collections for enrichment analysis.

``bulk.enrichment`` has always required a GMT. This module produces one from
GO, KEGG, and Reactome so a project does not have to supply its own, while
keeping the enrichment engine unchanged: the output is an ordinary GMT, so ORA
and GSEA both continue to run through clusterProfiler on one snapshot.

Fetching happens over the network, which no other TxSuite stage does. Isolating
it here is the point -- the collection becomes a workflow artifact with a dated
provenance record, so downstream stages stay offline and a later re-run can be
explained rather than merely repeated.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from txsuite.runtime import TxSuiteError

GENESET_SOURCES = ("go", "kegg", "reactome")
GENESET_KEYTYPES = ("symbol", "ensembl", "entrez")
GENESET_SPECIES = ("human", "mouse", "rat", "fly", "worm", "yeast", "zebrafish")
GO_ASPECTS = (
    "biological_process",
    "molecular_function",
    "cellular_component",
    "all",
)


def genesets_command(
    *,
    image: str,
    outdir: Path,
    sources: tuple[str, ...] = ("go",),
    species: str = "human",
    keytype: str = "symbol",
    aspect: str = "biological_process",
    min_size: int = 10,
    max_size: int = 500,
    min_mapped_fraction: float = 0.5,
) -> list[str]:
    """Build the container command that materializes a GMT."""

    if not sources:
        raise TxSuiteError("At least one gene-set source is required")
    unknown = sorted(set(sources) - set(GENESET_SOURCES))
    if unknown:
        raise TxSuiteError(
            f"Unknown gene-set source(s): {', '.join(unknown)}; "
            f"supported: {', '.join(GENESET_SOURCES)}"
        )
    if len(set(sources)) != len(sources):
        raise TxSuiteError("Gene-set sources must be unique")
    if species not in GENESET_SPECIES:
        raise TxSuiteError(f"Species must be one of: {', '.join(GENESET_SPECIES)}")
    if keytype not in GENESET_KEYTYPES:
        raise TxSuiteError(f"Keytype must be one of: {', '.join(GENESET_KEYTYPES)}")
    if aspect not in GO_ASPECTS:
        raise TxSuiteError(f"GO aspect must be one of: {', '.join(GO_ASPECTS)}")
    if min_size < 1 or max_size < min_size:
        raise TxSuiteError("Gene-set sizes must satisfy 1 <= min_size <= max_size")
    if not 0 < min_mapped_fraction <= 1:
        raise TxSuiteError("Minimum mapped fraction must be in (0, 1]")
    if not image.strip():
        raise TxSuiteError("Gene-set image cannot be empty")

    return [
        "docker",
        "run",
        "--rm",
        "--mount",
        f"type=bind,source={outdir.resolve()},target=/output",
        image,
        "python",
        "/opt/txsuite/fetch_genesets.py",
        "--outdir",
        "/output",
        "--sources",
        ",".join(sources),
        "--species",
        species,
        "--keytype",
        keytype,
        "--aspect",
        aspect,
        "--min-size",
        str(min_size),
        "--max-size",
        str(max_size),
        "--min-mapped-fraction",
        str(min_mapped_fraction),
    ]


def build_genesets_image(tag: str, *, run_dir: Path) -> None:
    """Build the gene-set fetching image from the packaged recipe."""

    import re
    import tempfile
    from importlib import resources

    from txsuite.runtime import run_command

    if not tag.strip():
        raise TxSuiteError("Image tag cannot be empty")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._:/-]*", tag):
        raise TxSuiteError("Image tag contains unsupported characters")
    package = resources.files("txsuite.resources.biodbs")
    with tempfile.TemporaryDirectory(prefix="txsuite-genesets-") as directory:
        context = Path(directory)
        for name in ("Dockerfile", "fetch_genesets.py"):
            (context / name).write_text(
                package.joinpath(name).read_text(encoding="utf-8"), encoding="utf-8"
            )
        run_command(
            ["docker", "build", "--tag", tag, str(context)],
            run_dir=run_dir,
            task="env.build.genesets",
            backend="docker",
            inputs={
                "base": (
                    "python:3.12-slim-bookworm@"
                    "sha256:d50fb7611f86d04a3b0471b46d7557818d88983fc3136726336b2a4c657aa30b"
                ),
                "biodbs": "0.4.1",
            },
            outputs={"image": tag},
            artifacts=[{"kind": "container-image", "label": "genesets", "path": tag}],
        )


__all__ = [
    "GENESET_KEYTYPES",
    "GENESET_SOURCES",
    "GENESET_SPECIES",
    "GO_ASPECTS",
    "build_genesets_image",
    "genesets_command",
]
