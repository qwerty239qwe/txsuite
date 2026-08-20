"""GATK RNA short variant discovery through pinned nf-core/rnavar.

The pipeline itself cannot run in CI, so these tests cover the parts TxSuite
owns: command construction, the validation that makes bad reference and
recalibration combinations fail before a run starts, and artifact resolution
against a synthetic copy of the rnavar output layout.
"""

from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG, ConfigError, load_config
from txsuite.project.config import load_project_config, parse_project_config
from txsuite.project.executor import ProjectExecutor
from txsuite.project.planner import plan_workflow
from txsuite.project.presets import scaffold_project_preset
from txsuite.project.provenance import RunBundle
from txsuite.project.results import StageState
from txsuite.project.adapters.bulk import bulk_rnavar_command
from txsuite.project.adapters.nfcore import (
    RNAVAR_1_3_0,
    get_nfcore_adapter,
    resolve_nfcore_artifacts,
)
from txsuite.project.registry import get_stage_spec, validate_stage_parameters
from txsuite.runtime import TxSuiteError
from txsuite.variants import rnavar_workflow_command, validate_rnavar_samplesheet


def _config():
    return load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))


class RnavarCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = _config()
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.samplesheet = self.root / "samplesheet.csv"
        self.samplesheet.write_text(
            "sample,fastq_1,fastq_2\nsample1,r1.fastq.gz,r2.fastq.gz\n", encoding="utf-8"
        )
        self.fasta = self.root / "genome.fa"
        self.fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        self.gtf = self.root / "genes.gtf"
        self.gtf.write_text("#gtf\n", encoding="utf-8")
        self.dbsnp = self.root / "dbsnp.vcf.gz"
        self.dbsnp.write_bytes(b"\x1f\x8b")
        self.outdir = self.root / "results"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _command(self, **kwargs):
        arguments = {
            "samplesheet": self.samplesheet,
            "outdir": self.outdir,
            "fasta": self.fasta,
            "gtf": self.gtf,
            "dbsnp": self.dbsnp,
        }
        arguments.update(kwargs)
        return rnavar_workflow_command(self.config, **arguments)

    def test_command_pins_the_release_and_passes_the_reference(self) -> None:
        command = self._command()

        self.assertEqual(command[:3], ["nextflow", "run", "nf-core/rnavar"])
        self.assertEqual(command[command.index("-r") + 1], "1.3.0")
        self.assertEqual(command[command.index("-profile") + 1], "docker")
        self.assertEqual(
            command[command.index("--input") + 1], str(self.samplesheet.resolve())
        )
        self.assertEqual(command[command.index("--fasta") + 1], str(self.fasta.resolve()))
        self.assertEqual(command[command.index("--dbsnp") + 1], str(self.dbsnp.resolve()))
        self.assertNotIn("--genome", command)
        self.assertNotIn("--skip_baserecalibration", command)
        self.assertNotIn("--tools", command)

    def test_igenomes_key_replaces_the_explicit_reference(self) -> None:
        command = self._command(genome="GRCh38", fasta=None, gtf=None)
        self.assertEqual(command[command.index("--genome") + 1], "GRCh38")
        self.assertNotIn("--fasta", command)

    def test_reference_selection_is_exclusive_and_complete(self) -> None:
        for kwargs in (
            {"fasta": None, "gtf": None},
            {"genome": "GRCh38"},
            {"gtf": None},
            {"fasta": None},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(TxSuiteError):
                self._command(**kwargs)

    def test_recalibration_requires_known_sites(self) -> None:
        with self.assertRaisesRegex(TxSuiteError, "requires --dbsnp"):
            self._command(dbsnp=None)

        skipped = self._command(dbsnp=None, skip_baserecalibration=True)
        self.assertIn("--skip_baserecalibration", skipped)

    def test_dbsnp_survives_skipped_recalibration_but_known_indels_do_not(self) -> None:
        # rnavar hands dbSNP to HaplotypeCaller for rsID annotation, which has
        # nothing to do with BQSR, so skipping recalibration must not reject it.
        command = self._command(skip_baserecalibration=True)
        self.assertIn("--skip_baserecalibration", command)
        self.assertEqual(command[command.index("--dbsnp") + 1], str(self.dbsnp.resolve()))

        indels = self.root / "known_indels.vcf.gz"
        indels.write_bytes(b"\x1f\x8b")
        with self.assertRaisesRegex(TxSuiteError, "Known indels are unused"):
            self._command(
                dbsnp=None, known_indels=indels, skip_baserecalibration=True
            )

    def test_unknown_tool_error_names_the_offending_value(self) -> None:
        with self.assertRaisesRegex(TxSuiteError, "bcftools"):
            self._command(tools=("bcftools",))

    def test_annotation_tools_require_their_caches(self) -> None:
        with self.assertRaisesRegex(TxSuiteError, "snpEff annotation requires"):
            self._command(tools=("snpeff",))
        with self.assertRaisesRegex(TxSuiteError, "VEP annotation requires"):
            self._command(tools=("vep",))
        with self.assertRaises(TxSuiteError):
            self._command(tools=("bcftools",))
        with self.assertRaises(TxSuiteError):
            self._command(tools=("vep", "vep"))

        cache = self.root / "snpeff_cache"
        cache.mkdir()
        command = self._command(tools=("snpeff",), snpeff_cache=cache)
        self.assertEqual(command[command.index("--tools") + 1], "snpeff")
        self.assertEqual(
            command[command.index("--snpeff_cache") + 1], str(cache.resolve())
        )

    def test_prebuilt_index_artifacts_are_forwarded(self) -> None:
        star_index = self.root / "star_index"
        star_index.mkdir()
        fai = self.root / "genome.fa.fai"
        fai.write_text("chr1\t4\t6\t60\t61\n", encoding="utf-8")
        sequence_dictionary = self.root / "genome.dict"
        sequence_dictionary.write_text("@HD\tVN:1.6\n", encoding="utf-8")

        command = self._command(
            star_index=star_index,
            fasta_fai=fai,
            sequence_dictionary=sequence_dictionary,
        )
        self.assertEqual(
            command[command.index("--star_index") + 1], str(star_index.resolve())
        )
        self.assertEqual(command[command.index("--fasta_fai") + 1], str(fai.resolve()))
        self.assertEqual(
            command[command.index("--dict") + 1], str(sequence_dictionary.resolve())
        )

    def test_gvcf_flag_is_opt_in(self) -> None:
        self.assertNotIn("--generate_gvcf", self._command())
        self.assertIn("--generate_gvcf", self._command(generate_gvcf=True))


class RnavarSamplesheetTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _sheet(self, text: str) -> Path:
        path = self.root / "samplesheet.csv"
        path.write_text(text, encoding="utf-8")
        return path

    def test_fastq_rows_are_accepted(self) -> None:
        path = self._sheet(
            "sample,fastq_1,fastq_2\na,a1.fq.gz,a2.fq.gz\nb,b1.fq.gz,b2.fq.gz\n"
        )
        self.assertEqual(validate_rnavar_samplesheet(path), 2)

    def test_other_entry_points_are_rejected_explicitly(self) -> None:
        path = self._sheet("sample,fastq_1,bam\na,a1.fq.gz,a.bam\n")
        with self.assertRaisesRegex(TxSuiteError, "only the FASTQ entry point"):
            validate_rnavar_samplesheet(path)

    def test_structural_problems_are_rejected(self) -> None:
        for text in (
            "sample\na\n",
            "sample,fastq_1\n",
            "sample,fastq_1\na,a1.fq.gz\na,a2.fq.gz\n",
            "sample,fastq_1\n,a1.fq.gz\n",
            "sample,fastq_1\nbad name,a1.fq.gz\n",
        ):
            with self.subTest(text=text), self.assertRaises(TxSuiteError):
                validate_rnavar_samplesheet(self._sheet(text))


class RnavarAdapterTests(unittest.TestCase):
    """Resolve artifacts against a synthetic copy of the rnavar output layout."""

    def _results(self, root: Path, samples: tuple[str, ...]) -> Path:
        results = (root / "results").resolve()
        for sample in samples:
            sample_dir = results / "variant_calling" / sample
            sample_dir.mkdir(parents=True)
            (sample_dir / f"{sample}.haplotypecaller.filtered.vcf.gz").write_bytes(b"\x1f\x8b")
            (sample_dir / f"{sample}.haplotypecaller.filtered.vcf.gz.tbi").write_bytes(b"\x00")
        multiqc = results / "reports" / "multiqc"
        multiqc.mkdir(parents=True)
        (multiqc / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
        return results

    def test_variant_directory_resolves_for_any_sample_count(self) -> None:
        for samples in (("sample1",), ("sample1", "sample2", "sample3")):
            with self.subTest(samples=len(samples)), tempfile.TemporaryDirectory() as d:
                results = self._results(Path(d), samples)

                resolved = resolve_nfcore_artifacts(
                    results,
                    adapter=RNAVAR_1_3_0,
                    artifact_types=(
                        "bulk.rnavar-results",
                        "bulk.variant-calls",
                        "qc.multiqc-report",
                    ),
                )

                self.assertEqual(
                    resolved["bulk.variant-calls"].path, results / "variant_calling"
                )
                self.assertEqual(
                    resolved["qc.multiqc-report"].path,
                    results / "reports" / "multiqc" / "multiqc_report.html",
                )
                self.assertEqual(
                    resolved["bulk.variant-calls"].evidence["release"], "1.3.0"
                )

    def test_missing_variant_calling_directory_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            results = Path(directory) / "results"
            results.mkdir()
            with self.assertRaises(TxSuiteError):
                resolve_nfcore_artifacts(
                    results,
                    adapter=RNAVAR_1_3_0,
                    artifact_types=("bulk.variant-calls",),
                )

    def test_adapter_is_registered_and_release_pinned(self) -> None:
        adapter = get_nfcore_adapter(pipeline="nf-core/rnavar", release="1.3.0")
        self.assertIs(adapter, RNAVAR_1_3_0)
        with self.assertRaises(TxSuiteError):
            get_nfcore_adapter(pipeline="nf-core/rnavar", release="1.2.0")


class RnavarStageTests(unittest.TestCase):
    def test_stage_declares_optional_reference_inputs(self) -> None:
        spec = get_stage_spec("bulk.rnavar")

        self.assertEqual(spec.modality, "bulk")
        self.assertEqual(set(spec.inputs), {"samplesheet"})
        self.assertEqual(
            set(spec.optional_inputs), {"star_index", "fasta_fai", "dict"}
        )
        self.assertEqual(
            spec.optional_inputs["star_index"],
            get_stage_spec("bulk.star-reference").outputs["star_index"],
        )
        self.assertEqual(spec.outputs["variants"], "bulk.variant-calls")
        self.assertEqual(spec.output_policies["variants"].kind, "directory")
        self.assertTrue(spec.supports_resume)

    def test_annotation_tool_parameter_is_validated(self) -> None:
        spec = get_stage_spec("bulk.rnavar")
        defaults = validate_stage_parameters(spec, {})
        self.assertEqual(defaults["tools"], ())
        self.assertFalse(defaults["generate_gvcf"])

        self.assertEqual(
            validate_stage_parameters(spec, {"tools": ["vep"]})["tools"], ("vep",)
        )
        for value in (["bcftools"], ["vep", "vep"], "vep"):
            with self.subTest(value=value), self.assertRaises(TxSuiteError):
                validate_stage_parameters(spec, {"tools": value})

    def test_adapter_matches_the_library_helper(self) -> None:
        config = _config()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            samplesheet = root / "samplesheet.csv"
            samplesheet.write_text(
                "sample,fastq_1,fastq_2\nsample1,r1.fq.gz,r2.fq.gz\n", encoding="utf-8"
            )
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
            gtf = root / "genes.gtf"
            gtf.write_text("#gtf\n", encoding="utf-8")
            outdir = root / "results"

            context = {
                "global_config": config,
                "resolved_inputs": {"samplesheet": samplesheet},
                "params": {
                    "fasta": fasta,
                    "gtf": gtf,
                    "skip_baserecalibration": True,
                },
                "outdir": outdir,
            }
            self.assertEqual(
                bulk_rnavar_command(context),
                rnavar_workflow_command(
                    config,
                    samplesheet=samplesheet,
                    outdir=outdir,
                    fasta=fasta,
                    gtf=gtf,
                    skip_baserecalibration=True,
                ),
            )


if __name__ == "__main__":
    unittest.main()


class RnavarPlanningTests(unittest.TestCase):
    """Planning must pin the right pipeline and hand postflight the right adapter.

    Both were driven by hardcoded stage tables that listed only the two original
    launchers, so a new launcher silently planned with no postflight metadata and
    its deferred artifacts would never have been resolved after a run.
    """

    def _plan(self, profile: str = "docker"):
        target = Path(tempfile.mkdtemp()) / "project"
        scaffold_project_preset("bulk-rnavar", target)
        workflow = load_project_config(target / "workflow.toml")
        if profile != "docker":
            workflow = parse_project_config(
                {
                    "schema_version": 1,
                    "project": {
                        "id": "bulk_rnavar_example",
                        "modality": "bulk",
                        "output_root": "results",
                    },
                    "execution": {"profile": profile, "resume": True},
                    "workflow": {
                        "stages": [
                            {
                                "id": "variants",
                                "uses": "bulk.rnavar",
                                "inputs": {"samplesheet": "samplesheet.csv"},
                                "params": {
                                    "fasta": "genome.fa",
                                    "gtf": "genes.gtf",
                                    "skip_baserecalibration": True,
                                },
                            }
                        ]
                    },
                },
                source_path=target / "workflow.toml",
            )
        return plan_workflow(workflow, DEFAULT_CONFIG)

    def test_postflight_names_the_pinned_rnavar_release(self) -> None:
        variants = self._plan().stage("variants")

        self.assertEqual(variants.postflight["adapter"], "nfcore")
        self.assertEqual(variants.postflight["pipeline"], "nf-core/rnavar")
        self.assertEqual(variants.postflight["release"], "1.3.0")
        self.assertEqual(
            variants.postflight["artifacts"]["variants"], "bulk.variant-calls"
        )
        self.assertEqual(
            dict(variants.postflight["policies"]["variants"]),
            {"kind": "directory", "non_empty": True},
        )

    def test_pipeline_pins_are_recorded_for_provenance(self) -> None:
        variants = self._plan().stage("variants")
        self.assertEqual(variants.pins["pipelines.variants.name"], "nf-core/rnavar")
        self.assertEqual(variants.pins["pipelines.variants.release"], "1.3.0")

    def test_native_stages_declare_no_postflight(self) -> None:
        reference = self._plan().stage("ref")
        self.assertEqual(dict(reference.postflight), {})
        self.assertIsNotNone(reference.outputs["star_index"].path)

    def test_nextflow_stages_plan_under_apptainer(self) -> None:
        variants = self._plan(profile="apptainer").stage("variants")
        self.assertEqual(variants.command[variants.command.index("-profile") + 1], "apptainer")


class RnavarExecutorTests(unittest.TestCase):
    """Postflight must resolve the variant directory and materialize downstream argv."""

    def _bundle(self, root: Path, plan: object) -> RunBundle:
        return RunBundle.create(
            root / "runs",
            run_id="rnavar-run",
            workflow={"name": "rnavar"},
            resolved_config={"profile": "test"},
            command_plan=plan,
        )

    def test_variant_directory_is_resolved_for_a_multi_sample_run(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workspace = root / "work"
            workspace.mkdir()
            results_root = workspace / "rnavar-results"
            downstream_output = workspace / "downstream.txt"

            raw_script = (
                "from pathlib import Path; import sys; "
                "root=Path(sys.argv[1]); "
                "[(root/'variant_calling'/s).mkdir(parents=True, exist_ok=True) "
                "for s in ('sample1', 'sample2')]; "
                "[(root/'variant_calling'/s/(s+'.haplotypecaller.filtered.vcf.gz'))"
                ".write_text('vcf') for s in ('sample1', 'sample2')]"
            )
            plan = {
                "stages": [
                    {
                        "id": "variants",
                        "command": [sys.executable, "-c", raw_script, str(results_root)],
                        "cwd": str(workspace),
                        "required_outputs": [],
                        "outputs": {"variants": "${variants.variants}"},
                        "output_artifacts": {
                            "variants": {
                                "type": "bulk.variant-calls",
                                "path": None,
                                "reference": "${variants.variants}",
                                "explicit": False,
                            }
                        },
                        "postflight": {
                            "adapter": "nfcore",
                            "pipeline": "nf-core/rnavar",
                            "release": "1.3.0",
                            "results_root": str(results_root),
                            "artifacts": {"variants": "bulk.variant-calls"},
                            "policies": {
                                "variants": {"kind": "directory", "non_empty": True}
                            },
                            "overrides": {},
                        },
                        "command_state": "resolved",
                    },
                    {
                        "id": "report",
                        "command": [
                            sys.executable,
                            "-c",
                            "from pathlib import Path; import sys; "
                            "source=Path(sys.argv[1].split('=', 1)[1]); "
                            "Path(sys.argv[2]).write_text("
                            "'\\n'.join(sorted(p.name for p in source.iterdir())))",
                            "input=${variants.variants}",
                            str(downstream_output),
                        ],
                        "cwd": str(workspace),
                        "required_outputs": [str(downstream_output)],
                        "inputs": {"variants": "${variants.variants}"},
                        "outputs": {"listing": str(downstream_output)},
                        "output_artifacts": {
                            "listing": {
                                "type": "test.listing",
                                "path": str(downstream_output),
                                "explicit": True,
                            }
                        },
                        "postflight": {},
                        "command_state": "deferred",
                    },
                ]
            }

            result = ProjectExecutor(self._bundle(root, plan)).execute()

            self.assertTrue(result.succeeded)
            # The downstream stage received the directory, not an unresolved token.
            self.assertEqual(
                downstream_output.read_text(encoding="utf-8").splitlines(),
                ["sample1", "sample2"],
            )
            artifact = result.stages[0].artifacts["variants"]
            self.assertEqual(
                artifact["path"], str((results_root / "variant_calling").resolve())
            )
            self.assertEqual(artifact["evidence"]["adapter"], "nfcore-rnavar-1.3.0")
            self.assertEqual(artifact["evidence"]["release"], "1.3.0")
            report_command = result.stages[1].attempts[0].command
            self.assertFalse(
                any("${variants.variants}" in item for item in report_command)
            )


class VariantConfigPinTests(unittest.TestCase):
    def test_blank_rnavar_pin_or_star_image_is_rejected(self) -> None:
        # Packaged defaults always backfill a missing key, so the failure mode a
        # project file can actually produce is a blank override, not an absent one.
        for section, key, blanked in (
            ("images", "star", {"star": ""}),
            ("pipelines", "variants", {"variants": {"name": "", "release": "1.3.0"}}),
        ):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as d:
                config = copy.deepcopy(DEFAULT_CONFIG)
                config[section].update(blanked)
                path = Path(d) / "txsuite.toml"
                path.write_text(_toml(config), encoding="utf-8")
                with self.assertRaises(ConfigError):
                    load_config(path, user_path=Path(d) / "missing.toml")

    def test_defaults_pin_rnavar_and_the_star_image(self) -> None:
        config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        self.assertEqual(config["pipelines"]["variants"]["name"], "nf-core/rnavar")
        self.assertEqual(config["pipelines"]["variants"]["release"], "1.3.0")
        self.assertTrue(config["images"]["star"].startswith("txsuite/star:"))


def _toml(config: dict) -> str:
    lines = ["[execution]", f'profile = "{config["execution"]["profile"]}"', "", "[images]"]
    for name, value in config["images"].items():
        lines.append(f'{name} = "{value}"')
    for name, pipeline in config["pipelines"].items():
        lines.extend(
            ["", f"[pipelines.{name}]", f'name = "{pipeline["name"]}"',
             f'release = "{pipeline["release"]}"']
        )
    return "\n".join(lines) + "\n"
