from __future__ import annotations

import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path

from txsuite.bulk import enrichment_command
from txsuite.config import DEFAULT_CONFIG
from txsuite.project.registry import (
    OutputPolicy,
    StageSpec,
    get_stage_spec,
    list_stage_specs,
    validate_stage_parameters,
)
from txsuite.runtime import TxSuiteError
from txsuite.single_cell import (
    analysis_command,
    pseudobulk_command,
    pseudobulk_workflow_command,
)


class StageRegistryTests(unittest.TestCase):
    def test_scanpy_registry_forwards_analysis_controls(self):
        spec = get_stage_spec("single-cell.scanpy")
        params = validate_stage_parameters(spec, {
            "counts_layer": "raw", "stop_after": "pca", "skip_umap": True,
            "skip_markers": True, "n_pcs": 12, "n_neighbors": 8,
            "integration": "harmony", "batch_column": "batch", "doublets": "score",
        })
        with tempfile.TemporaryDirectory() as temporary:
            command = spec.build_command({
                "global_config": DEFAULT_CONFIG, "resolved_inputs": {"input": Path(temporary) / "input.h5ad"},
                "params": params, "outdir": Path(temporary) / "out", "check_inputs": False,
            })
        for flag, value in (("--counts-layer", "raw"), ("--stop-after", "pca"),
                            ("--n-pcs", "12"), ("--integration", "harmony"), ("--doublets", "score")):
            self.assertEqual(command[command.index(flag) + 1], value)
        self.assertIn("--skip-umap", command)
        self.assertIn("--skip-markers", command)

    def test_registry_contains_the_wave_one_stages_in_stable_order(self) -> None:
        expected = (
            "bulk.rnaseq",
            "bulk.de",
            "bulk.enrichment",
            "single-cell.scrnaseq",
            "single-cell.scanpy",
            "single-cell.pseudobulk",
            "single-cell.pseudobulk-de",
        )
        specs = list_stage_specs()
        self.assertEqual(tuple(spec.uses for spec in specs), expected)
        self.assertEqual(get_stage_spec("bulk.de").id, "bulk.de")
        self.assertEqual(
            tuple(spec.id for spec in list_stage_specs("bulk")), expected[:3]
        )
        with self.assertRaises(FrozenInstanceError):
            specs[0].maturity = "ready"  # type: ignore[misc]
        for spec in specs:
            self.assertEqual(set(spec.output_policies), set(spec.outputs))
            for policy in spec.output_policies.values():
                self.assertIsInstance(policy, OutputPolicy)
                self.assertIn(policy.kind, {"file", "directory"})
            for artifact_type in spec.outputs.values():
                self.assertNotRegex(artifact_type, r"[*?\[\]/\\]")
        self.assertTrue(get_stage_spec("single-cell.scanpy").output_policies["h5ad"].non_empty)
        self.assertEqual(
            get_stage_spec("single-cell.scrnaseq").output_policies["matrix"].kind,
            "file",
        )
        self.assertTrue(
            get_stage_spec("single-cell.scrnaseq")
            .output_policies["matrix"]
            .non_empty
        )
        with self.assertRaises(TypeError):
            specs[0].output_policies["results"] = OutputPolicy("file")  # type: ignore[index]
        with self.assertRaises(FrozenInstanceError):
            specs[0].output_policies["results"].non_empty = False  # type: ignore[misc]
        self.assertEqual(
            get_stage_spec("bulk.de").outputs["de_results"],
            get_stage_spec("bulk.enrichment").inputs["de_results"],
        )
        self.assertEqual(
            get_stage_spec("single-cell.pseudobulk").outputs["counts"],
            get_stage_spec("bulk.de").inputs["counts"],
        )
        self.assertIn("multiqc_report", get_stage_spec("bulk.rnaseq").outputs)
        self.assertIn(
            "multiqc_report", get_stage_spec("single-cell.scrnaseq").outputs
        )

    def test_parameters_are_defaulted_normalized_and_strict(self) -> None:
        de = get_stage_spec("bulk.de")
        params = validate_stage_parameters(
            de,
            {
                "design": "condition",
                "reference": "control",
                "test": "treated",
                "covariates": ["batch"],
            },
        )
        self.assertEqual(params["method"], "deseq2")
        self.assertEqual(params["covariates"], ("batch",))
        with self.assertRaisesRegex(TxSuiteError, "Unknown parameter"):
            validate_stage_parameters(de, {**params, "typo": True})
        with self.assertRaisesRegex(TxSuiteError, "Missing required parameter"):
            validate_stage_parameters(de, {})
        with self.assertRaisesRegex(TxSuiteError, "reference and test must differ"):
            validate_stage_parameters(
                de,
                {"design": "condition", "reference": "same", "test": "same"},
            )

    def test_stage_defaults_are_recursively_immutable_copies(self) -> None:
        supplied = {"nested": {"items": [{"threshold": 1}]}}
        spec = StageSpec(
            id="test.deep-freeze",
            modality="bulk",
            maturity="ready",
            inputs={},
            outputs={},
            defaults=supplied,
            validators={},
            command_factory=lambda context: [],
        )
        supplied["nested"]["items"][0]["threshold"] = 2
        nested = spec.defaults["nested"]
        self.assertEqual(nested["items"][0]["threshold"], 1)
        self.assertIsInstance(nested["items"], tuple)
        with self.assertRaises(TypeError):
            nested["items"][0]["threshold"] = 3

    def test_output_policy_defaults_remain_backward_compatible_and_partial_is_rejected(self) -> None:
        spec = StageSpec(
            id="test.output-policy-default",
            modality="bulk",
            maturity="ready",
            inputs={},
            outputs={"results": "test.results", "table": "test.table"},
            defaults={},
            validators={},
            command_factory=lambda context: [],
        )
        self.assertEqual(spec.output_policies["results"].kind, "directory")
        self.assertEqual(spec.output_policies["table"].kind, "file")
        positional = StageSpec(
            "test.output-policy-positional",
            "bulk",
            "ready",
            {},
            {"results": "test.results"},
            {},
            {},
            lambda context: [],
            ("tool",),
            ("images.test",),
            True,
            "test.output-policy-positional",
        )
        self.assertEqual(positional.required_executables, ("tool",))
        self.assertEqual(positional.required_images, ("images.test",))
        self.assertTrue(positional.supports_resume)
        with self.assertRaisesRegex(ValueError, "exactly match"):
            StageSpec(
                id="test.output-policy-partial",
                modality="bulk",
                maturity="ready",
                inputs={},
                outputs={"one": "test.one", "two": "test.two"},
                output_policies={"one": OutputPolicy("file")},
                defaults={},
                validators={},
                command_factory=lambda context: [],
            )

    def test_command_factory_uses_existing_bulk_builder(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            samplesheet = root / "samples.csv"
            samplesheet.write_text(
                "sample,fastq_1,fastq_2,strandedness\n"
                "a,a_R1.fastq.gz,a_R2.fastq.gz,auto\n",
                encoding="utf-8",
            )
            command = get_stage_spec("bulk.rnaseq").build_command(
                {
                    "global_config": DEFAULT_CONFIG,
                    "resolved_inputs": {"samplesheet": samplesheet},
                    "params": {},
                    "outdir": root / "results",
                    "resume": True,
                }
            )
            self.assertEqual(command[:5], ["nextflow", "run", "nf-core/rnaseq", "-r", "3.26.0"])
            self.assertEqual(command[-1], "-resume")

    def test_unknown_stage_fails_explicitly(self) -> None:
        with self.assertRaisesRegex(TxSuiteError, "Unknown project stage"):
            get_stage_spec("bulk.unknown")

    def test_command_builders_can_defer_future_input_checks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            future_de = root / "future" / "de.tsv"
            future_genesets = root / "future" / "sets.gmt"
            future_matrix = root / "future" / "matrix"
            future_h5ad = root / "future" / "data.h5ad"

            enrichment = enrichment_command(
                image="txsuite/bulk-r:test",
                de_results=future_de,
                genesets=future_genesets,
                mode="ora",
                outdir=root / "enrichment",
                check_inputs=False,
            )
            scanpy = analysis_command(
                image="txsuite/single-cell-python:test",
                input_path=future_matrix,
                outdir=root / "scanpy",
                min_genes=1,
                min_cells=1,
                max_mito_pct=20,
                resolution=1.0,
                check_inputs=False,
            )
            aggregate = pseudobulk_command(
                image="txsuite/single-cell-python:test",
                h5ad=future_h5ad,
                outdir=root / "pseudobulk",
                sample_column="sample",
                design="condition",
                check_inputs=False,
            )
            workflow = pseudobulk_workflow_command(
                DEFAULT_CONFIG,
                h5ad=future_h5ad,
                outdir=root / "pseudobulk-de",
                sample_column="sample",
                design="condition",
                reference="control",
                test="treated",
                check_inputs=False,
            )

            self.assertEqual(enrichment[:3], ["docker", "run", "--rm"])
            self.assertEqual(scanpy[:3], ["docker", "run", "--rm"])
            self.assertEqual(aggregate[:3], ["docker", "run", "--rm"])
            self.assertEqual(workflow[:2], ["nextflow", "run"])
            with self.assertRaisesRegex(TxSuiteError, "Enrichment mode"):
                enrichment_command(
                    image="txsuite/bulk-r:test",
                    de_results=future_de,
                    genesets=future_genesets,
                    mode="invalid",
                    outdir=root / "invalid",
                    check_inputs=False,
                )

    def test_command_builders_remain_strict_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            missing = root / "missing"
            calls = (
                lambda: enrichment_command(
                    image="txsuite/bulk-r:test",
                    de_results=missing / "de.tsv",
                    genesets=missing / "sets.gmt",
                    mode="ora",
                    outdir=root / "enrichment",
                ),
                lambda: analysis_command(
                    image="txsuite/single-cell-python:test",
                    input_path=missing / "matrix",
                    outdir=root / "scanpy",
                    min_genes=1,
                    min_cells=1,
                    max_mito_pct=20,
                    resolution=1.0,
                ),
                lambda: pseudobulk_command(
                    image="txsuite/single-cell-python:test",
                    h5ad=missing / "data.h5ad",
                    outdir=root / "pseudobulk",
                    sample_column="sample",
                    design="condition",
                ),
                lambda: pseudobulk_workflow_command(
                    DEFAULT_CONFIG,
                    h5ad=missing / "data.h5ad",
                    outdir=root / "pseudobulk-de",
                    sample_column="sample",
                    design="condition",
                    reference="control",
                    test="treated",
                ),
            )
            for call in calls:
                with self.subTest(call=call), self.assertRaisesRegex(
                    TxSuiteError, "does not exist"
                ):
                    call()

    def test_stage_context_passes_deferred_input_policy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            command = get_stage_spec("single-cell.pseudobulk").build_command(
                {
                    "global_config": DEFAULT_CONFIG,
                    "resolved_inputs": {"h5ad": root / "future.h5ad"},
                    "params": {"sample_column": "sample", "design": "condition"},
                    "outdir": root / "out",
                    "resume": False,
                    "check_inputs": False,
                }
            )
            self.assertIn("/opt/txsuite/single_cell.py", command)


if __name__ == "__main__":
    unittest.main()
