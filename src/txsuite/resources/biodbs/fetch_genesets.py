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


def translate_symbols(
    symbols: list[str], *, keytype: str, species: str
) -> dict[str, str]:
    """Map fetched gene identifiers into ``keytype``. Networked."""

    from biodbs import biomart_convert_ids

    attribute = _BIOMART_ATTRIBUTE[keytype]
    dataset = _BIOMART_DATASET[species]
    try:
        converted = biomart_convert_ids(
            ids=symbols,
            from_type="external_gene_name",
            to_type=attribute,
            dataset=dataset,
        )
    except Exception as exc:  # noqa: BLE001
        raise GeneSetError(f"gene identifier translation failed: {exc}") from exc
    frame = converted.to_pandas() if hasattr(converted, "to_pandas") else converted
    mapping: dict[str, str] = {}
    for row in frame.itertuples(index=False):
        source_value, target_value = str(row[0]).strip(), str(row[1]).strip()
        if source_value and target_value and target_value.lower() != "nan":
            mapping.setdefault(source_value, target_value)
    return mapping


def apply_translation(
    collections: dict[str, dict], mapping: dict[str, str]
) -> tuple[dict[str, dict], dict[str, int]]:
    """Rewrite gene members through ``mapping``, counting what was lost.

    Unmapped members are dropped rather than passed through: leaving an
    untranslated identifier in a collection keyed by another namespace would
    silently never match the differential-expression table.
    """

    translated: dict[str, dict] = {}
    seen: set[str] = set()
    mapped: set[str] = set()
    for key, record in collections.items():
        members = []
        for gene in record["genes"]:
            seen.add(gene)
            target = mapping.get(gene)
            if target is None:
                continue
            mapped.add(gene)
            members.append(target)
        translated[key] = {**record, "genes": sorted(set(members))}
    return translated, {"input": len(seen), "mapped": len(mapped)}


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
    counts: dict[str, int],
    versions: dict[str, str],
    parameters: dict[str, object],
    biodbs_version: str,
    min_mapped_fraction: float,
    fetched_at: str,
) -> dict[str, Path]:
    """Write the GMT, the mapping report, and the provenance record."""

    total, mapped = counts.get("input", 0), counts.get("mapped", 0)
    fraction = 1.0 if total == 0 else mapped / total
    if total and fraction < min_mapped_fraction:
        raise GeneSetError(
            f"only {mapped}/{total} gene identifiers ({fraction:.1%}) mapped to "
            f"{parameters['keytype']}; below the {min_mapped_fraction:.0%} floor. "
            "Check that keytype matches the differential-expression table"
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
        writer.writerow(["metric", "value"])
        writer.writerow(["keytype", parameters["keytype"]])
        writer.writerow(["input_identifiers", total])
        writer.writerow(["mapped_identifiers", mapped])
        writer.writerow(["unmapped_identifiers", total - mapped])
        writer.writerow(["mapped_fraction", f"{fraction:.4f}"])
        for source in sorted(per_source):
            writer.writerow([f"terms_{source}", per_source[source]])

    written["provenance"].write_text(
        json.dumps(
            {
                "fetched_at": fetched_at,
                "biodbs_version": biodbs_version,
                "sources": versions,
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
        counts = {"input": 0, "mapped": 0}
        if arguments.keytype != "symbol":
            symbols = sorted({g for r in collections.values() for g in r["genes"]})
            mapping = translate_symbols(
                symbols, keytype=arguments.keytype, species=arguments.species
            )
            collections, counts = apply_translation(collections, mapping)
        else:
            identifiers = {g for r in collections.values() for g in r["genes"]}
            counts = {"input": len(identifiers), "mapped": len(identifiers)}

        import biodbs

        write_outputs(
            arguments.outdir,
            collections,
            counts=counts,
            versions=versions,
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
