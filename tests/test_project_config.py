from __future__ import annotations

import json
import re
import tempfile
import unittest
from copy import deepcopy
from dataclasses import FrozenInstanceError
from pathlib import Path

from txsuite.project import (
    ArtifactReference,
    WorkflowConfigError,
    load_project_config,
    load_workflow_schema,
    parse_artifact_reference,
    parse_project_config,
    redact_sensitive,
)

try:
    import jsonschema
except ImportError:  # Optional test dependency; TxSuite itself stays dependency-free.
    jsonschema = None


VALID_TOML = """\
schema_version = 1

[project]
id = "example_project"
modality = "bulk"
output_root = "results"

[execution]
profile = "docker"
resume = true

[[workflow.stages]]
id = "align"
uses = "bulk.rnaseq"

[workflow.stages.inputs]
samplesheet = "data/samples.csv"
reads = { paired = ["data/R1.fastq.gz", "data/R2.fastq.gz"] }

[workflow.stages.params]
threads = 4
reference_like_literal = "${literal.value}"
credentials = { password = "top-secret", label = "visible" }

[workflow.stages.outputs]
counts = "artifacts/counts.tsv"
bam = "artifacts/aligned.bam"

[[workflow.stages]]
id = "differential"
uses = "bulk.de"
depends_on = ["align"]

[workflow.stages.inputs]
counts = "${align.counts}"
nested = { primary = ["${align.bam}"] }

[workflow.stages.params]
design = "~ condition"

[workflow.stages.outputs]
report = "artifacts/report.html"
"""


def valid_document() -> dict[str, object]:
    return {
        "schema_version": 1,
        "project": {
            "id": "example",
            "modality": "bulk",
            "output_root": "results",
        },
        "execution": {"profile": "docker", "resume": False},
        "workflow": {
            "stages": [
                {
                    "id": "prepare",
                    "uses": "bulk.rnaseq",
                    "inputs": {"source": "data/input.h5ad"},
                    "params": {},
                    "outputs": {"matrix": "artifacts/matrix.h5ad"},
                }
            ]
        },
    }


class ProjectConfigTests(unittest.TestCase):
    def test_loads_and_resolves_canonical_toml(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = root / "config" / "project.toml"
            config_path.parent.mkdir()
            config_path.write_text(VALID_TOML, encoding="utf-8")

            workflow = load_project_config(config_path)

            self.assertEqual(workflow.schema_version, 1)
            self.assertEqual(workflow.source_path, config_path.resolve())
            self.assertEqual(workflow.source_dir, config_path.parent.resolve())
            self.assertEqual(workflow.project.id, "example_project")
            self.assertEqual(workflow.project.modality, "bulk")
            self.assertEqual(
                workflow.project.output_root,
                (config_path.parent / "results").resolve(),
            )
            self.assertEqual(workflow.execution.profile, "docker")
            self.assertTrue(workflow.execution.resume)
            self.assertEqual(
                workflow.stage("align").inputs["samplesheet"],
                (config_path.parent / "data/samples.csv").resolve(),
            )
            self.assertEqual(
                workflow.stage("align").inputs["reads"]["paired"][1],
                (config_path.parent / "data/R2.fastq.gz").resolve(),
            )
            self.assertEqual(
                workflow.stage("align").outputs["counts"],
                (config_path.parent / "artifacts/counts.tsv").resolve(),
            )

    def test_artifact_references_are_typed_recursively_but_params_are_literal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "project.toml"
            path.write_text(VALID_TOML, encoding="utf-8")
            workflow = load_project_config(path)

        stage = workflow.stage("differential")
        self.assertEqual(stage.inputs["counts"], ArtifactReference("align", "counts"))
        self.assertEqual(
            stage.inputs["nested"]["primary"][0],
            ArtifactReference("align", "bam"),
        )
        self.assertEqual(
            workflow.stage("align").params["reference_like_literal"],
            "${literal.value}",
        )

    def test_resolved_models_are_deeply_immutable(self) -> None:
        workflow = parse_project_config(
            valid_document(), source_path=Path("somewhere/project.toml")
        )
        with self.assertRaises(FrozenInstanceError):
            workflow.schema_version = 2  # type: ignore[misc]
        with self.assertRaises(TypeError):
            workflow.stage("prepare").params["new"] = True  # type: ignore[index]
        with self.assertRaises(AttributeError):
            workflow.stages.append(workflow.stages[0])  # type: ignore[attr-defined]

    def test_canonical_representation_is_json_serializable_and_deterministic(self) -> None:
        workflow = parse_project_config(
            valid_document(), source_path=Path("somewhere/project.toml")
        )
        encoded = workflow.to_json()
        decoded = json.loads(encoded)

        self.assertEqual(decoded, workflow.to_dict())
        self.assertEqual(encoded, workflow.to_json())
        self.assertEqual(decoded["schema_version"], 1)
        self.assertEqual(decoded["workflow"]["stages"][0]["depends_on"], [])
        self.assertIsInstance(
            decoded["workflow"]["stages"][0]["inputs"]["source"], str
        )

    def test_redacts_sensitive_keys_at_arbitrary_depth(self) -> None:
        value = {
            "token": "one",
            "nested": [
                {"API-Key": "two", "safe": "shown"},
                {"credentials": {"password": "three"}},
            ],
        }
        self.assertEqual(
            redact_sensitive(value),
            {
                "token": "***REDACTED***",
                "nested": [
                    {"API-Key": "***REDACTED***", "safe": "shown"},
                    {"credentials": "***REDACTED***"},
                ],
            },
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "project.toml"
            path.write_text(VALID_TOML, encoding="utf-8")
            redacted = load_project_config(path).to_dict(redact=True)
        credentials = redacted["workflow"]["stages"][0]["params"]["credentials"]
        self.assertEqual(credentials, "***REDACTED***")

    def test_rejects_unknown_keys_at_every_schema_layer(self) -> None:
        cases = (
            (("surprise",), "configuration"),
            (("project", "owner"), "project"),
            (("execution", "queue"), "execution"),
            (("workflow", "name"), "workflow"),
            (("workflow", "stages", 0, "command"), "workflow.stages[0]"),
        )
        for path, expected in cases:
            with self.subTest(path=path):
                document = valid_document()
                current: object = document
                for part in path[:-1]:
                    current = current[part]  # type: ignore[index]
                current[path[-1]] = "not allowed"  # type: ignore[index]
                with self.assertRaisesRegex(WorkflowConfigError, re.escape(expected)):
                    parse_project_config(
                        document, source_path=Path("project.toml")
                    )

    def test_rejects_missing_required_tables_and_fields(self) -> None:
        for key in ("schema_version", "project", "execution", "workflow"):
            with self.subTest(key=key):
                document = valid_document()
                del document[key]
                with self.assertRaisesRegex(WorkflowConfigError, "required"):
                    parse_project_config(document, source_path="project.toml")

        document = valid_document()
        del document["project"]["output_root"]  # type: ignore[index]
        with self.assertRaisesRegex(WorkflowConfigError, "project.output_root"):
            parse_project_config(document, source_path="project.toml")

    def test_rejects_unsupported_or_wrong_type_schema_versions(self) -> None:
        for value in (2, "1", True):
            with self.subTest(value=value):
                document = valid_document()
                document["schema_version"] = value
                with self.assertRaisesRegex(WorkflowConfigError, "schema_version"):
                    parse_project_config(document, source_path="project.toml")

    def test_rejects_unsupported_modality_and_profile(self) -> None:
        for modality in ("proteomics", "single_cell"):
            with self.subTest(modality=modality):
                document = valid_document()
                document["project"]["modality"] = modality  # type: ignore[index]
                with self.assertRaisesRegex(WorkflowConfigError, "project.modality"):
                    parse_project_config(document, source_path="project.toml")

        for profile in ("local", "slurm"):
            with self.subTest(profile=profile):
                document = valid_document()
                document["execution"]["profile"] = profile  # type: ignore[index]
                with self.assertRaisesRegex(WorkflowConfigError, "execution.profile"):
                    parse_project_config(document, source_path="project.toml")

    def test_rejects_duplicate_stage_and_dependency_ids(self) -> None:
        document = valid_document()
        first = document["workflow"]["stages"][0]  # type: ignore[index]
        document["workflow"]["stages"].append(dict(first))  # type: ignore[index]
        with self.assertRaisesRegex(WorkflowConfigError, "duplicate workflow stage ID"):
            parse_project_config(document, source_path="project.toml")

        document = valid_document()
        document["workflow"]["stages"].append(  # type: ignore[index]
            {
                "id": "consume",
                "uses": "bulk.de",
                "depends_on": ["prepare", "prepare"],
            }
        )
        with self.assertRaisesRegex(WorkflowConfigError, "duplicate stage IDs"):
            parse_project_config(document, source_path="project.toml")

    def test_rejects_malformed_artifact_references(self) -> None:
        for reference in ("${prepare}", "prefix-${prepare.matrix}", "${.matrix}"):
            with self.subTest(reference=reference):
                document = valid_document()
                document["workflow"]["stages"][0]["inputs"]["source"] = reference  # type: ignore[index]
                with self.assertRaisesRegex(WorkflowConfigError, "malformed artifact"):
                    parse_project_config(document, source_path="project.toml")

        self.assertIsNone(parse_artifact_reference("ordinary/path"))

    def test_rejects_references_to_unknown_stages(self) -> None:
        document = valid_document()
        document["workflow"]["stages"].append(  # type: ignore[index]
            {
                "id": "consume",
                "uses": "bulk.de",
                "inputs": {"value": "${missing.matrix}"},
            }
        )
        with self.assertRaisesRegex(WorkflowConfigError, "unknown stage"):
            parse_project_config(document, source_path="project.toml")

    def test_artifact_name_validation_is_deferred_and_reference_implies_dependency(self) -> None:
        document = valid_document()
        document["workflow"]["stages"] = [  # type: ignore[index]
            {"id": "raw", "uses": "bulk.rnaseq"},
            {
                "id": "consume",
                "uses": "bulk.de",
                "inputs": {"matrix": "${raw.matrix}"},
            },
        ]
        workflow = parse_project_config(document, source_path="project.toml")

        raw = workflow.stage("raw")
        consume = workflow.stage("consume")
        self.assertEqual(dict(raw.outputs), {})
        self.assertEqual(consume.depends_on, ())
        self.assertEqual(consume.dependencies, ("raw",))
        self.assertEqual(
            consume.artifact_references,
            (ArtifactReference("raw", "matrix"),),
        )

    def test_rejects_unknown_dependencies_self_dependencies_and_cycles(self) -> None:
        document = valid_document()
        document["workflow"]["stages"][0]["depends_on"] = ["missing"]  # type: ignore[index]
        with self.assertRaisesRegex(WorkflowConfigError, "unknown stage"):
            parse_project_config(document, source_path="project.toml")

        document = valid_document()
        document["workflow"]["stages"][0]["depends_on"] = ["prepare"]  # type: ignore[index]
        with self.assertRaisesRegex(WorkflowConfigError, "depend on itself"):
            parse_project_config(document, source_path="project.toml")

        document = valid_document()
        document["workflow"]["stages"][0]["depends_on"] = ["consume"]  # type: ignore[index]
        document["workflow"]["stages"].append(  # type: ignore[index]
            {
                "id": "consume",
                "uses": "bulk.de",
                "depends_on": ["prepare"],
            }
        )
        with self.assertRaisesRegex(WorkflowConfigError, "dependency cycle"):
            parse_project_config(document, source_path="project.toml")

        document = valid_document()
        document["workflow"]["stages"] = [  # type: ignore[index]
            {
                "id": "first",
                "uses": "bulk.rnaseq",
                "inputs": {"value": "${second.output}"},
            },
            {
                "id": "second",
                "uses": "bulk.de",
                "inputs": {"value": "${first.output}"},
            },
        ]
        with self.assertRaisesRegex(WorkflowConfigError, "dependency cycle"):
            parse_project_config(document, source_path="project.toml")

    def test_rejects_non_json_toml_params_and_non_finite_numbers(self) -> None:
        document = valid_document()
        document["workflow"]["stages"][0]["params"] = {"when": object()}  # type: ignore[index]
        with self.assertRaisesRegex(WorkflowConfigError, "JSON-serializable"):
            parse_project_config(document, source_path="project.toml")

        document = valid_document()
        document["workflow"]["stages"][0]["params"] = {"bad": float("inf")}  # type: ignore[index]
        with self.assertRaisesRegex(WorkflowConfigError, "finite number"):
            parse_project_config(document, source_path="project.toml")

    def test_wraps_missing_file_and_toml_decode_errors(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(WorkflowConfigError, "cannot read"):
                load_project_config(root / "missing.toml")

            malformed = root / "malformed.toml"
            malformed.write_text("schema_version = [", encoding="utf-8")
            with self.assertRaisesRegex(WorkflowConfigError, "cannot read"):
                load_project_config(malformed)

    def test_packaged_schema_locks_the_canonical_shape(self) -> None:
        schema = load_workflow_schema()
        self.assertEqual(schema["properties"]["schema_version"]["const"], 1)
        self.assertIn(
            "does not distinguish",
            schema["properties"]["schema_version"]["description"],
        )
        self.assertFalse(schema["additionalProperties"])
        self.assertEqual(
            set(schema["required"]),
            {"schema_version", "project", "execution", "workflow"},
        )
        stage = schema["$defs"]["stage"]
        self.assertFalse(stage["additionalProperties"])
        self.assertEqual(set(stage["required"]), {"id", "uses"})
        self.assertEqual(
            set(schema["$defs"]["execution"]["properties"]["profile"]["enum"]),
            {"docker", "apptainer"},
        )
        self.assertEqual(stage["properties"]["uses"]["pattern"], "\\S")

    @unittest.skipIf(jsonschema is None, "optional jsonschema package is not installed")
    def test_runtime_and_json_schema_have_table_driven_parity(self) -> None:
        schema = load_workflow_schema()
        validator = jsonschema.Draft202012Validator(schema)  # type: ignore[union-attr]

        referenced = valid_document()
        del referenced["workflow"]["stages"][0]["outputs"]  # type: ignore[index]
        referenced["workflow"]["stages"].append(  # type: ignore[index]
            {
                "id": "consume",
                "uses": "bulk.de",
                "inputs": {"matrix": "${prepare.matrix}"},
            }
        )

        local_profile = valid_document()
        local_profile["execution"]["profile"] = "local"  # type: ignore[index]
        blank_uses = valid_document()
        blank_uses["workflow"]["stages"][0]["uses"] = " \t"  # type: ignore[index]
        blank_root = valid_document()
        blank_root["project"]["output_root"] = "  "  # type: ignore[index]
        blank_output = valid_document()
        blank_output["workflow"]["stages"][0]["outputs"]["matrix"] = "\t"  # type: ignore[index]
        blank_input = valid_document()
        blank_input["workflow"]["stages"][0]["inputs"]["source"] = " "  # type: ignore[index]
        malformed_reference = valid_document()
        malformed_reference["workflow"]["stages"][0]["inputs"]["source"] = (  # type: ignore[index]
            "prefix-${prepare.matrix}"
        )
        unknown_key = valid_document()
        unknown_key["project"]["owner"] = "analyst"  # type: ignore[index]
        modality_alias = valid_document()
        modality_alias["project"]["modality"] = "single_cell"  # type: ignore[index]

        cases = (
            ("canonical", valid_document(), True),
            ("exact artifact reference", referenced, True),
            ("local profile", local_profile, False),
            ("whitespace uses", blank_uses, False),
            ("whitespace output root", blank_root, False),
            ("whitespace output", blank_output, False),
            ("whitespace input path", blank_input, False),
            ("malformed artifact reference", malformed_reference, False),
            ("unknown key", unknown_key, False),
            ("modality alias", modality_alias, False),
        )

        for label, document, expected in cases:
            with self.subTest(label=label):
                schema_valid = validator.is_valid(document)
                try:
                    parse_project_config(
                        deepcopy(document), source_path="project.toml"
                    )
                except WorkflowConfigError:
                    runtime_valid = False
                else:
                    runtime_valid = True
                self.assertEqual(schema_valid, expected)
                self.assertEqual(runtime_valid, expected)

    @unittest.skipIf(jsonschema is None, "optional jsonschema package is not installed")
    def test_schema_documents_json_number_toml_integer_limitation(self) -> None:
        document = valid_document()
        document["schema_version"] = 1.0
        schema = load_workflow_schema()
        validator = jsonschema.Draft202012Validator(schema)  # type: ignore[union-attr]

        self.assertTrue(validator.is_valid(document))
        with self.assertRaisesRegex(WorkflowConfigError, "schema_version"):
            parse_project_config(document, source_path="project.toml")


if __name__ == "__main__":
    unittest.main()
