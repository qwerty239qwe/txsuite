from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from txsuite.project.adapters.nfcore import (
    ArtifactRule,
    NfCoreAdapter,
    SCRNASEQ_4_2_0,
    resolve_nfcore_artifact,
    resolve_nfcore_artifacts,
)
from txsuite.runtime import TxSuiteError


class NfCoreArtifactTests(unittest.TestCase):
    def setUp(self) -> None:
        self.adapter = NfCoreAdapter(
            name="fixture-1.2.3",
            pipeline="nf-core/fixture",
            release="1.2.3",
            artifacts={
                "fixture.table": ArtifactRule(("tables/result.tsv", "alt/result.tsv")),
                "fixture.results": ArtifactRule((".",), "directory"),
            },
        )

    def test_unique_resolution_carries_adapter_and_release_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table = root / "tables" / "result.tsv"
            table.parent.mkdir()
            table.write_text("result\n", encoding="utf-8")
            resolved = resolve_nfcore_artifact(
                root, "fixture.table", adapter=self.adapter
            )
            self.assertEqual(resolved.path, table.resolve())
            self.assertEqual(resolved.evidence["adapter"], "fixture-1.2.3")
            self.assertEqual(resolved.evidence["release"], "1.2.3")
            self.assertEqual(resolved.evidence["source"], "adapter")

    def test_zero_and_ambiguous_matches_fail(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(TxSuiteError, "No 'fixture.table' artifact"):
                resolve_nfcore_artifact(root, "fixture.table", adapter=self.adapter)
            for parent in ("tables", "alt"):
                path = root / parent / "result.tsv"
                path.parent.mkdir(exist_ok=True)
                path.touch()
            with self.assertRaisesRegex(TxSuiteError, "Ambiguous 'fixture.table'"):
                resolve_nfcore_artifact(root, "fixture.table", adapter=self.adapter)

    def test_explicit_override_is_deterministic_and_recorded(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for parent in ("tables", "alt"):
                path = root / parent / "result.tsv"
                path.parent.mkdir(exist_ok=True)
                path.touch()
            resolved = resolve_nfcore_artifact(
                root,
                "fixture.table",
                adapter=self.adapter,
                override="alt/result.tsv",
            )
            self.assertEqual(resolved.path, (root / "alt/result.tsv").resolve())
            self.assertEqual(resolved.evidence["source"], "override")
            self.assertNotIn("pattern", resolved.evidence)

    def test_plural_resolution_and_release_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table = root / "chosen.tsv"
            table.touch()
            resolved = resolve_nfcore_artifacts(
                root,
                adapter=self.adapter,
                artifact_types=("fixture.table", "fixture.results"),
                overrides={"fixture.table": table},
            )
            self.assertEqual(set(resolved), {"fixture.table", "fixture.results"})
            with self.assertRaisesRegex(TxSuiteError, "not '9.9.9'"):
                resolve_nfcore_artifact(
                    root,
                    "fixture.table",
                    adapter=self.adapter,
                    release="9.9.9",
                    override=table,
                )

    def test_scrnaseq_resolves_only_the_combined_anndata_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            combined = root / "simpleaf" / "mtx_conversions" / "combined_matrix.h5ad"
            combined.parent.mkdir(parents=True)
            combined.touch()
            qcatch = root / "simpleaf" / "qcatch" / "sample.h5ad"
            qcatch.parent.mkdir()
            qcatch.touch()

            resolved = resolve_nfcore_artifact(
                root, "single-cell.matrix", adapter=SCRNASEQ_4_2_0
            )
            self.assertEqual(resolved.path, combined.resolve())

            second = root / "starsolo" / "mtx_conversions" / "combined_matrix.h5ad"
            second.parent.mkdir(parents=True)
            second.touch()
            with self.assertRaisesRegex(TxSuiteError, "Ambiguous 'single-cell.matrix'"):
                resolve_nfcore_artifact(
                    root, "single-cell.matrix", adapter=SCRNASEQ_4_2_0
                )

    def test_artifact_rules_reject_parent_traversal(self) -> None:
        with self.assertRaisesRegex(ValueError, "stay below"):
            ArtifactRule(("tables/../outside.tsv",))

    def test_overrides_must_remain_below_the_results_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            root = base / "results"
            root.mkdir()
            outside = base / "outside.tsv"
            outside.touch()
            for override in ("../outside.tsv", outside.resolve()):
                with self.subTest(override=override), self.assertRaisesRegex(
                    TxSuiteError, "cannot contain|escapes"
                ):
                    resolve_nfcore_artifact(
                        root,
                        "fixture.table",
                        adapter=self.adapter,
                        override=override,
                    )

    def test_symlink_artifacts_cannot_escape_the_results_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            root = base / "results"
            root.mkdir()
            outside = base / "outside"
            outside.mkdir()
            (outside / "result.tsv").touch()
            (root / "tables").symlink_to(outside, target_is_directory=True)

            with self.assertRaisesRegex(TxSuiteError, "escapes"):
                resolve_nfcore_artifact(root, "fixture.table", adapter=self.adapter)
            with self.assertRaisesRegex(TxSuiteError, "escapes"):
                resolve_nfcore_artifact(
                    root,
                    "fixture.table",
                    adapter=self.adapter,
                    override="tables/result.tsv",
                )


if __name__ == "__main__":
    unittest.main()
