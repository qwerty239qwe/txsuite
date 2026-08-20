"""Materialize gene-set collections into a GMT plus an audit trail.

This is the only TxSuite step that reaches the network at run time. It exists so
that everything downstream does not have to: the gene sets are fetched once,
written to a GMT that becomes an ordinary workflow artifact, and accompanied by
a provenance record. A result produced six months later can then be explained,
because the collection it used is on disk and the fetch is dated.

Fetching is confined to :func:`fetch_collections`. Everything else -- ID
translation bookkeeping, GMT rendering, the mapping report -- is pure, so it is
tested without a network.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

SOURCES = ("go", "kegg", "reactome")
KEYTYPES = ("symbol", "ensembl", "entrez")
SPECIES = ("human", "mouse", "rat", "fly", "worm", "yeast", "zebrafish")
GO_ASPECTS = ("biological_process", "molecular_function", "cellular_component", "all")

# biomaRt attribute names for each keytype we expose, plus the dataset per
# species. Kept as data so an unsupported combination fails on a lookup here
# rather than inside a query.
_BIOMART_ATTRIBUTE = {
    "symbol": "external_gene_name",
    "ensembl": "ensembl_gene_id",
    "entrez": "entrezgene_id",
}
# Each source returns members in its own namespace -- biodbs documents KEGG as
# Entrez, GO (QuickGO) as UniProt, and Reactome as gene symbols. Translating
# every source as if it were symbols silently produces a collection that never
# matches the differential-expression table.
_SOURCE_ID_TYPE = {
    "go": "uniprot_gn_id",
    "kegg": "entrezgene_id",
    "reactome": "external_gene_name",
}
_BIOMART_DATASET = {
    "human": "hsapiens_gene_ensembl",
    "mouse": "mmusculus_gene_ensembl",
    "rat": "rnorvegicus_gene_ensembl",
    "fly": "dmelanogaster_gene_ensembl",
    "worm": "celegans_gene_ensembl",
    "yeast": "scerevisiae_gene_ensembl",
    "zebrafish": "drerio_gene_ensembl",
}


class GeneSetError(RuntimeError):
    """A gene-set collection could not be fetched, translated, or written."""


def fetch_collections(
    sources: list[str],
    *,
    species: str,
    aspect: str,
    min_size: int,
    max_size: int,
) -> tuple[dict[str, dict], dict[str, str]]:
    """Fetch each source, returning ``{term_id: record}`` and source versions.

    The only networked function in this module.
    """

    from biodbs.analysis import fetch_gmt  # imported here so the rest stays pure

    collections: dict[str, dict] = {}
    versions: dict[str, str] = {}
    for source in sources:
        try:
            fetched = fetch_gmt(
                name=f"txsuite-{source}",
                database=source,
                species=species,
                aspect=aspect,
                min_term_size=min_size,
                max_term_size=max_size,
            )
        except Exception as exc:  # noqa: BLE001 - upstream raises many API types
            raise GeneSetError(f"could not fetch {source} gene sets: {exc}") from exc
        if not fetched:
            raise GeneSetError(f"{source} returned no gene sets for species {species!r}")
        for term_id, pathway in fetched.items():
            key = f"{source}:{term_id}"
            collections[key] = {
                "source": source,
                "term_id": term_id,
                "name": getattr(pathway, "name", term_id) or term_id,
                "genes": sorted(getattr(pathway, "genes", ()) or ()),
            }
        versions[source] = str(getattr(next(iter(fetched.values())), "database", source))
    return collections, versions


def translate_members(
    identifiers: list[str], *, from_type: str, to_type: str, species: str
) -> dict[str, str]:
    """Map identifiers from one namespace to another. Networked."""

    from biodbs import biomart_convert_ids

    dataset = _BIOMART_DATASET[species]
    try:
        converted = biomart_convert_ids(
            ids=identifiers, from_type=from_type, to_type=to_type, dataset=dataset
        )
    except Exception as exc:  # noqa: BLE001
        raise GeneSetError(
            f"translation from {from_type} to {to_type} failed: {exc}"
        ) from exc
    frame = converted.to_pandas() if hasattr(converted, "to_pandas") else converted
    mapping: dict[str, str] = {}
    for row in frame.itertuples(index=False):
        source_value, target_value = str(row[0]).strip(), str(row[1]).strip()
        if source_value and target_value and target_value.lower() != "nan":
            mapping.setdefault(source_value, target_value)
    return mapping


def apply_translation(
    collections: dict[str, dict], mappings: dict[str, dict[str, str]]
) -> tuple[dict[str, dict], dict[str, dict[str, int]]]:
    """Rewrite members through each source's own mapping, counting the losses.

    ``mappings`` is keyed by source, because the namespace a source returns is a
    property of that source. An unmapped member is dropped rather than passed
    through: keeping an identifier from the wrong namespace would look like a
    gene set member and never match anything.
    """

    translated: dict[str, dict] = {}
    seen: dict[str, set[str]] = {}
    mapped: dict[str, set[str]] = {}
    for key, record in collections.items():
        source = record["source"]
        mapping = mappings.get(source)
        seen.setdefault(source, set())
        mapped.setdefault(source, set())
        members = []
        for gene in record["genes"]:
            seen[source].add(gene)
            if mapping is None:  # already in the requested namespace
                mapped[source].add(gene)
                members.append(gene)
                continue
            target = mapping.get(gene)
            if target is None:
                continue
            mapped[source].add(gene)
            members.append(target)
        translated[key] = {**record, "genes": sorted(set(members))}
    counts = {
        source: {"input": len(seen[source]), "mapped": len(mapped[source])}
        for source in seen
    }
    return translated, counts


def render_gmt(collections: dict[str, dict]) -> str:
    """Render a GMT deterministically.

    Terms and members are sorted so the artifact hashes identically for the same
    content; the run bundle's resume logic compares hashes, and an unordered
    dump would look changed on every fetch.
    """

    lines = []
    for key in sorted(collections):
        record = collections[key]
        if not record["genes"]:
            continue
        description = f"{record['source']}|{record['name']}"
        lines.append("\t".join([record["term_id"], description, *sorted(record["genes"])]))
    if not lines:
        raise GeneSetError("no gene sets survived filtering and translation")
    return "\n".join(lines) + "\n"


def write_outputs(
    outdir: Path,
    collections: dict[str, dict],
    *,
    counts: dict[str, dict[str, int]],
    sources_used: dict[str, str],
    parameters: dict[str, object],
    biodbs_version: str,
    min_mapped_fraction: float,
    fetched_at: str,
) -> dict[str, Path]:
    """Write the GMT, the per-source mapping report, and the provenance record."""

    total = sum(entry["input"] for entry in counts.values())
    mapped = sum(entry["mapped"] for entry in counts.values())
    fraction = 1.0 if total == 0 else mapped / total
    # Checked per source as well as overall: one broken namespace among three
    # can stay above a global floor while contributing nothing.
    for source in sorted(counts):
        entry = counts[source]
        if not entry["input"]:
            continue
        source_fraction = entry["mapped"] / entry["input"]
        if source_fraction < min_mapped_fraction:
            raise GeneSetError(
                f"{source}: only {entry['mapped']}/{entry['input']} identifiers "
                f"({source_fraction:.1%}) mapped to {parameters['keytype']}, below "
                f"the {min_mapped_fraction:.0%} floor. {source} returns "
                f"{_SOURCE_ID_TYPE.get(source, 'unknown')} identifiers; check that "
                "keytype matches the differential-expression table"
            )

    outdir.mkdir(parents=True, exist_ok=True)
    written = {
        "gmt": outdir / "genesets.gmt",
        "mapping": outdir / "id-mapping.tsv",
        "provenance": outdir / "sources.json",
    }
    written["gmt"].write_text(render_gmt(collections), encoding="utf-8")

    per_source: dict[str, int] = {}
    for record in collections.values():
        per_source[record["source"]] = per_source.get(record["source"], 0) + 1
    with written["mapping"].open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["source", "native_id_type", "terms", "input", "mapped", "unmapped", "fraction"]
        )
        for source in sorted(counts):
            entry = counts[source]
            unmapped = entry["input"] - entry["mapped"]
            ratio = 1.0 if not entry["input"] else entry["mapped"] / entry["input"]
            writer.writerow(
                [
                    source,
                    _SOURCE_ID_TYPE.get(source, "unknown"),
                    per_source.get(source, 0),
                    entry["input"],
                    entry["mapped"],
                    unmapped,
                    f"{ratio:.4f}",
                ]
            )

    written["provenance"].write_text(
        json.dumps(
            {
                "fetched_at": fetched_at,
                "biodbs_version": biodbs_version,
                # The APIs do not expose a release identifier, so this records
                # which source produced each collection, not a pinned version.
                "sources": sources_used,
                "source_id_types": {s: _SOURCE_ID_TYPE.get(s, "unknown") for s in counts},
                "parameters": parameters,
                "terms": sum(per_source.values()),
                "mapped_fraction": round(fraction, 4),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--sources", required=True, help="comma-separated: go,kegg,reactome")
    parser.add_argument("--species", default="human", choices=SPECIES)
    parser.add_argument("--keytype", default="symbol", choices=KEYTYPES)
    parser.add_argument("--aspect", default="biological_process", choices=GO_ASPECTS)
    parser.add_argument("--min-size", type=int, default=10)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--min-mapped-fraction", type=float, default=0.5)
    arguments = parser.parse_args(argv)

    sources = [item.strip() for item in arguments.sources.split(",") if item.strip()]
    unknown = sorted(set(sources) - set(SOURCES))
    if unknown or not sources:
        print(
            f"error: sources must be a non-empty subset of {', '.join(SOURCES)}",
            file=sys.stderr,
        )
        return 1

    try:
        collections, versions = fetch_collections(
            sources,
            species=arguments.species,
            aspect=arguments.aspect,
            min_size=arguments.min_size,
            max_size=arguments.max_size,
        )
        target = _BIOMART_ATTRIBUTE[arguments.keytype]
        mappings: dict[str, dict[str, str]] = {}
        for source in sources:
            native = _SOURCE_ID_TYPE[source]
            if native == target:
                continue  # already in the requested namespace
            members = sorted(
                {
                    gene
                    for record in collections.values()
                    if record["source"] == source
                    for gene in record["genes"]
                }
            )
            if members:
                mappings[source] = translate_members(
                    members,
                    from_type=native,
                    to_type=target,
                    species=arguments.species,
                )
        collections, counts = apply_translation(collections, mappings)

        import biodbs

        write_outputs(
            arguments.outdir,
            collections,
            counts=counts,
            sources_used=versions,
            parameters={
                "sources": sources,
                "species": arguments.species,
                "keytype": arguments.keytype,
                "aspect": arguments.aspect,
                "min_size": arguments.min_size,
                "max_size": arguments.max_size,
            },
            biodbs_version=getattr(biodbs, "__version__", "unknown"),
            min_mapped_fraction=arguments.min_mapped_fraction,
            fetched_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )
    except GeneSetError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
