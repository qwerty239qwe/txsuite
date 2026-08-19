"""Materialized gene-set collections.

This is the only stage that reaches the network, so the tests are deliberately
split: the fetch itself is one isolated function that is never called here, and
everything that shapes the artifacts -- translation bookkeeping, GMT rendering,
the mapping floor, provenance -- is pure and exercised directly.
"""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG, load_config
from txsuite.genesets import (
    GENESET_SOURCES,
    build_genesets_image,
    genesets_command,
)
from txsuite.project.adapters.bulk import bulk_genesets_command
from txsuite.project.config import parse_project_config
from txsuite.project.planner import plan_workflow
from txsuite.project.registry import get_stage_spec, validate_stage_parameters
from txsuite.runtime import TxSuiteError


def _script():
    path = Path(
        str(resources.files("txsuite.resources.biodbs").joinpath("fetch_genesets.py"))
    )
    spec = importlib.util.spec_from_file_location("txsuite_fetch_genesets", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _collection(source: str, term: str, name: str, genes: list[str]) -> dict:
    return {f"{source}:{term}": {"source": source, "term_id": term, "name": name, "genes": sorted(genes)}}


class GenesetsCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.directory = tempfile.TemporaryDirectory()
        self.outdir = Path(self.directory.name) / "genesets"
        self.outdir.mkdir()

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _command(self, **kwargs):
        arguments = {"image": self.config["images"]["genesets"], "outdir": self.outdir}
        arguments.update(kwargs)
        return genesets_command(**arguments)

    def test_sources_and_settings_reach_the_container(self) -> None:
        command = self._command(
            sources=("go", "kegg"), species="mouse", keytype="ensembl"
        )
        self.assertEqual(command[command.index("--sources") + 1], "go,kegg")
        self.assertEqual(command[command.index("--species") + 1], "mouse")
        self.assertEqual(command[command.index("--keytype") + 1], "ensembl")
        self.assertIn("/opt/txsuite/fetch_genesets.py", command)
        self.assertEqual(command[command.index("--outdir") + 1], "/output")

    def test_invalid_settings_are_rejected(self) -> None:
        for kwargs in (
            {"sources": ()},
            {"sources": ("go", "go")},
            {"sources": ("msigdb",)},
            {"species": "axolotl"},
            {"keytype": "refseq"},
            {"aspect": "everything"},
            {"min_size": 0},
            {"min_size": 100, "max_size": 10},
            {"min_mapped_fraction": 0.0},
            {"min_mapped_fraction": 1.5},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                self._command(**kwargs)

    def test_adapter_matches_the_library_helper(self) -> None:
        context = {
            "global_config": self.config,
            "resolved_inputs": {},
            "params": {"sources": ["kegg", "reactome"], "keytype": "entrez"},
            "outdir": self.outdir,
        }
        self.assertEqual(
            bulk_genesets_command(context),
            self._command(sources=("kegg", "reactome"), keytype="entrez"),
        )

    def test_image_build_rejects_bad_tags_before_docker_runs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            run_dir = Path(directory) / "run"
            for tag in ("", "  ", "bad tag"):
                with self.subTest(tag=tag), self.assertRaises(TxSuiteError):
                    build_genesets_image(tag, run_dir=run_dir)
            self.assertFalse(run_dir.exists())


class GenesetsStageTests(unittest.TestCase):
    def test_gmt_output_type_matches_what_enrichment_consumes(self) -> None:
        spec = get_stage_spec("bulk.genesets")
        # The whole design rests on this: the emitted GMT is the same artifact
        # type bulk.enrichment already takes, so ORA and GSEA are unchanged.
        self.assertEqual(
            spec.outputs["gmt"], get_stage_spec("bulk.enrichment").inputs["genesets"]
        )
        self.assertEqual(dict(spec.inputs), {})
        self.assertEqual(spec.required_images, ("images.genesets",))

    def test_defaults_and_validation(self) -> None:
        spec = get_stage_spec("bulk.genesets")
        defaults = validate_stage_parameters(spec, {})
        self.assertEqual(defaults["sources"], ("go",))
        self.assertEqual(defaults["species"], "human")
        self.assertEqual(defaults["keytype"], "symbol")
        self.assertEqual(
            validate_stage_parameters(spec, {"sources": ["kegg", "go"]})["sources"],
            ("kegg", "go"),
        )
        for params in ({"sources": []}, {"sources": ["go", "go"]}, {"sources": "go"},
                       {"sources": ["msigdb"]}, {"species": "axolotl"}):
            with self.subTest(params=params), self.assertRaises(TxSuiteError):
                validate_stage_parameters(spec, params)

    def test_registry_and_library_agree_on_sources(self) -> None:
        spec = get_stage_spec("bulk.genesets")
        for source in GENESET_SOURCES:
            with self.subTest(source=source):
                self.assertEqual(
                    validate_stage_parameters(spec, {"sources": [source]})["sources"],
                    (source,),
                )

    def test_plan_wires_the_gmt_into_enrichment(self) -> None:
        workflow = parse_project_config(
            {
                "schema_version": 1,
                "project": {"id": "g", "modality": "bulk", "output_root": "results"},
                "execution": {"profile": "docker", "resume": True},
                "workflow": {
                    "stages": [
                        {
                            "id": "genesets",
                            "uses": "bulk.genesets",
                            "params": {"sources": ["go"]},
                        },
                        {
                            "id": "enrichment",
                            "uses": "bulk.enrichment",
                            "inputs": {
                                "de_results": "de.tsv",
                                "genesets": "${genesets.gmt}",
                            },
                        },
                    ]
                },
            },
            source_path=Path("/work/workflow.toml"),
        )
        plan = plan_workflow(workflow, DEFAULT_CONFIG)
        genesets, enrichment = plan.stage("genesets"), plan.stage("enrichment")
        self.assertEqual(genesets.outputs["gmt"].path, genesets.outdir / "genesets.gmt")
        self.assertIn("genesets", enrichment.depends_on)
        self.assertIn(str(genesets.outputs["gmt"].path), " ".join(enrichment.command))


class GenesetArtifactTests(unittest.TestCase):
    """The pure half of the fetch script: rendering, mapping, provenance."""

    def setUp(self) -> None:
        self.module = _script()
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _write(self, collections, *, counts=None, keytype="symbol", floor=0.5):
        return self.module.write_outputs(
            self.root / "out",
            collections,
            counts=counts or {"go": {"input": 3, "mapped": 3}},
            sources_used={"go": "go"},
            parameters={"sources": ["go"], "keytype": keytype},
            biodbs_version="0.4.1",
            min_mapped_fraction=floor,
            fetched_at="2026-08-19T00:00:00+00:00",
        )

    def test_gmt_is_rendered_deterministically(self) -> None:
        # Sorting matters beyond tidiness: the run bundle skips a stage by
        # comparing artifact hashes, so unordered output would look changed on
        # every fetch and defeat resume.
        unsorted = {
            "go:GO:2": {"source": "go", "term_id": "GO:2", "name": "b", "genes": ["Z", "A"]},
            "go:GO:1": {"source": "go", "term_id": "GO:1", "name": "a", "genes": ["M", "B"]},
        }
        rendered = self.module.render_gmt(unsorted)
        self.assertEqual(
            rendered.splitlines(),
            ["GO:1\tgo|a\tB\tM", "GO:2\tgo|b\tA\tZ"],
        )
        shuffled = {k: unsorted[k] for k in reversed(list(unsorted))}
        self.assertEqual(self.module.render_gmt(shuffled), rendered)

    def test_empty_terms_are_dropped_and_an_empty_collection_fails(self) -> None:
        mixed = {
            **_collection("go", "GO:1", "kept", ["A"]),
            **_collection("go", "GO:2", "empty", []),
        }
        self.assertEqual(len(self.module.render_gmt(mixed).splitlines()), 1)
        with self.assertRaises(self.module.GeneSetError):
            self.module.render_gmt({**_collection("go", "GO:2", "empty", [])})

    def test_each_source_is_translated_from_its_own_namespace(self) -> None:
        # biodbs returns KEGG as Entrez, GO as UniProt and Reactome as symbols.
        # Translating all three as symbols silently yields sets that match
        # nothing, so the mapping is keyed by source.
        collections = {
            **_collection("kegg", "hsa1", "one", ["7157", "672"]),
            **_collection("go", "GO:1", "two", ["P04637"]),
            **_collection("reactome", "R-HSA-1", "three", ["TP53"]),
        }
        translated, counts = self.module.apply_translation(
            collections,
            {
                "kegg": {"7157": "ENSG00000141510"},
                "go": {"P04637": "ENSG00000141510"},
                # reactome absent: already in the requested namespace
            },
        )
        self.assertEqual(translated["kegg:hsa1"]["genes"], ["ENSG00000141510"])
        self.assertEqual(translated["go:GO:1"]["genes"], ["ENSG00000141510"])
        self.assertEqual(translated["reactome:R-HSA-1"]["genes"], ["TP53"])
        self.assertEqual(counts["kegg"], {"input": 2, "mapped": 1})
        self.assertEqual(counts["go"], {"input": 1, "mapped": 1})
        self.assertEqual(counts["reactome"], {"input": 1, "mapped": 1})

    def test_source_native_id_types_match_what_biodbs_returns(self) -> None:
        self.assertEqual(
            self.module._SOURCE_ID_TYPE,
            {
                "go": "uniprot_gn_id",
                "kegg": "entrezgene_id",
                "reactome": "external_gene_name",
            },
        )

    def test_one_broken_source_fails_even_when_the_others_carry_the_average(self) -> None:
        # A global floor alone would pass this: 210/310 overall is 68%, while
        # KEGG contributed nothing at all.
        with self.assertRaisesRegex(self.module.GeneSetError, "kegg: only 0/100"):
            self._write(
                _collection("kegg", "hsa1", "one", ["A"]),
                counts={
                    "kegg": {"input": 100, "mapped": 0},
                    "go": {"input": 210, "mapped": 210},
                },
                keytype="ensembl",
            )

    def test_a_poor_mapping_rate_fails_instead_of_enriching_on_a_fragment(self) -> None:
        collections = _collection("kegg", "hsa1", "one", ["A"])
        with self.assertRaisesRegex(self.module.GeneSetError, "below the"):
            self._write(
                collections,
                counts={"kegg": {"input": 100, "mapped": 20}},
                keytype="entrez",
            )
        written = self._write(
            collections, counts={"kegg": {"input": 100, "mapped": 80}}, keytype="entrez"
        )
        self.assertTrue(written["gmt"].is_file())

    def test_provenance_records_what_would_explain_a_changed_result(self) -> None:
        written = self._write(_collection("go", "GO:1", "a", ["A", "B", "C"]))
        record = json.loads(written["provenance"].read_text(encoding="utf-8"))
        self.assertEqual(record["fetched_at"], "2026-08-19T00:00:00+00:00")
        self.assertEqual(record["biodbs_version"], "0.4.1")
        self.assertEqual(record["sources"], {"go": "go"})
        self.assertEqual(record["terms"], 1)
        self.assertEqual(record["parameters"]["keytype"], "symbol")

    def test_mapping_report_records_the_loss_per_source(self) -> None:
        written = self._write(
            _collection("kegg", "hsa1", "one", ["A"]),
            counts={"kegg": {"input": 200, "mapped": 150}},
            keytype="entrez",
        )
        rows = [
            line.split("\t")
            for line in written["mapping"].read_text(encoding="utf-8").splitlines()
        ]
        self.assertEqual(
            rows[0],
            ["source", "native_id_type", "terms", "input", "mapped", "unmapped", "fraction"],
        )
        self.assertEqual(rows[1], ["kegg", "entrezgene_id", "1", "200", "150", "50", "0.7500"])

    def test_the_script_does_not_import_biodbs_at_module_scope(self) -> None:
        # Keeps the pure half testable without the dependency, and keeps the
        # network confined to the one function that needs it.
        source = (
            resources.files("txsuite.resources.biodbs")
            .joinpath("fetch_genesets.py")
            .read_text(encoding="utf-8")
        )
        for line in source.splitlines():
            self.assertFalse(
                line.startswith("import biodbs") or line.startswith("from biodbs"),
                f"module-scope biodbs import: {line}",
            )


if __name__ == "__main__":
    unittest.main()
