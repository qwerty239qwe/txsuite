"""Contrast expansion on ``bulk.de``.

A multi-group study should not need one hand-written stage per comparison. The
``contrasts`` parameter derives the comparison set from the design column's
levels while leaving the primary contrast, and every existing artifact, exactly
as a single-contrast run produces them.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.bulk import CONTRAST_MODES, differential_expression_command
from txsuite.config import load_config
from txsuite.project.adapters.bulk import bulk_de_command
from txsuite.project.registry import get_stage_spec, validate_stage_parameters
from txsuite.runtime import TxSuiteError


class ContrastCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.counts = self.root / "counts.tsv"
        self.counts.write_text("gene_id\ts1\ts2\ngene1\t5\t7\n", encoding="utf-8")
        self.metadata = self.root / "metadata.tsv"
        self.metadata.write_text("sample\tcondition\ns1\tcontrol\ns2\ttreated\n", encoding="utf-8")
        self.outdir = self.root / "de"

    def tearDown(self) -> None:
        self.directory.cleanup()

    def _command(self, method: str = "deseq2", **kwargs) -> list[str]:
        return differential_expression_command(
            method=method,
            image="txsuite/bulk-r:test",
            counts=self.counts,
            metadata=self.metadata,
            outdir=self.outdir,
            design="condition",
            reference="control",
            test="treated",
            **kwargs,
        )

    def test_mode_is_the_last_argument_and_defaults_to_single(self) -> None:
        self.assertEqual(self._command()[-1], "single")
        self.assertEqual(self._command(contrasts="all-pairs")[-1], "all-pairs")
        self.assertEqual(self._command(contrasts="vs-reference")[-1], "vs-reference")

    def test_every_method_forwards_the_mode(self) -> None:
        for method in ("deseq2", "edger", "limma"):
            with self.subTest(method=method):
                command = self._command(method=method, contrasts="all-pairs")
                self.assertEqual(command[-1], "all-pairs")

    def test_unknown_modes_are_rejected(self) -> None:
        with self.assertRaises(TxSuiteError):
            self._command(contrasts="everything")

    def test_formula_mode_cannot_expand_contrasts(self) -> None:
        with self.assertRaises(TxSuiteError):
            differential_expression_command(
                method="deseq2",
                image="txsuite/bulk-r:test",
                counts=self.counts,
                metadata=self.metadata,
                outdir=self.outdir,
                formula="~ condition",
                coefficient="conditiontreated",
                contrasts="all-pairs",
            )

    def test_adapter_forwards_the_registry_parameter(self) -> None:
        config = load_config(Path("does-not-exist.toml"), user_path=Path("missing.toml"))
        context = {
            "global_config": config,
            "resolved_inputs": {"counts": self.counts, "metadata": self.metadata},
            "params": {
                "design": "condition",
                "reference": "control",
                "test": "treated",
                "method": "deseq2",
                "contrasts": "vs-reference",
            },
            "outdir": self.outdir,
        }
        self.assertEqual(bulk_de_command(context)[-1], "vs-reference")


class ContrastRegistryTests(unittest.TestCase):
    def test_default_preserves_single_contrast_behaviour(self) -> None:
        spec = get_stage_spec("bulk.de")
        defaults = validate_stage_parameters(
            spec, {"design": "condition", "reference": "control", "test": "treated"}
        )
        self.assertEqual(defaults["contrasts"], "single")

    def test_every_documented_mode_validates(self) -> None:
        spec = get_stage_spec("bulk.de")
        for mode in CONTRAST_MODES:
            with self.subTest(mode=mode):
                params = validate_stage_parameters(
                    spec,
                    {
                        "design": "condition",
                        "reference": "control",
                        "test": "treated",
                        "contrasts": mode,
                    },
                )
                self.assertEqual(params["contrasts"], mode)
        with self.assertRaises(TxSuiteError):
            validate_stage_parameters(
                spec,
                {
                    "design": "condition",
                    "reference": "control",
                    "test": "treated",
                    "contrasts": "pairwise",
                },
            )

    def test_contrast_index_is_declared_and_verified(self) -> None:
        spec = get_stage_spec("bulk.de")
        self.assertEqual(
            spec.outputs["contrast_index"], "bulk.differential-expression-index"
        )
        policy = spec.output_policies["contrast_index"]
        self.assertEqual(policy.kind, "file")
        self.assertTrue(policy.non_empty)
        # The enrichment edge still consumes the primary contrast unchanged.
        self.assertEqual(
            spec.outputs["de_results"],
            get_stage_spec("bulk.enrichment").inputs["de_results"],
        )


class ContrastScriptTests(unittest.TestCase):
    """The R scripts are not runnable here, so assert their contract textually."""

    def _script(self, name: str) -> str:
        return resources.files("txsuite.resources.bulk_r").joinpath(name).read_text(
            encoding="utf-8"
        )

    def test_both_scripts_accept_and_validate_the_mode(self) -> None:
        for name, lengths in (
            ("deseq2.R", "c(6, 10, 12, 13)"),
            ("alternative_de.R", "c(11, 13, 14)"),
        ):
            with self.subTest(script=name):
                script = self._script(name)
                self.assertIn(lengths, script)
                self.assertIn(
                    "contrasts must be 'single', 'vs-reference', or 'all-pairs'", script
                )
                self.assertIn(
                    "formula mode has no design levels to expand", script
                )

    def test_both_scripts_write_the_index_and_per_contrast_tables(self) -> None:
        for name in ("deseq2.R", "alternative_de.R"):
            with self.subTest(script=name):
                script = self._script(name)
                self.assertIn('file.path(outdir, "contrasts.tsv")', script)
                self.assertIn('paste0("DE_", pair[[1]], "_vs_", pair[[2]], ".tsv")', script)
                self.assertIn("contrast_id", script)
                # The expansion is capped so a high-cardinality column cannot
                # silently launch a quadratic number of comparisons.
                self.assertIn("> 50L", script)

    def test_alternative_script_orders_levels_reference_first(self) -> None:
        # deseq2.R gets this from relevel(); alternative_de.R has no full-factor
        # model, so it must order the levels itself before expanding. Using raw
        # alphabetical order there once produced an inverted duplicate of the
        # primary contrast.
        script = self._script("alternative_de.R")
        self.assertIn(
            "design_levels <- c(reference, setdiff(observed_levels, reference))",
            script,
        )
        self.assertIn("relevel(metadata[[design]], ref = reference)", self._script("deseq2.R"))

    def test_python_and_r_agree_on_the_accepted_modes(self) -> None:
        for name in ("deseq2.R", "alternative_de.R"):
            with self.subTest(script=name):
                script = self._script(name)
                match = re.search(
                    r'contrast_mode %in% c\((.*?)\)', script, re.DOTALL
                )
                self.assertIsNotNone(match)
                modes = tuple(re.findall(r'"([^"]+)"', match.group(1)))
                self.assertEqual(modes, CONTRAST_MODES)


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is not installed")
class ContrastSetTests(unittest.TestCase):
    """Exercise the real contrast-set logic.

    ``expanded_contrasts`` is plain R with no Bioconductor dependency, so it can
    be extracted from the shipped script and evaluated directly. That keeps the
    comparison set — the part where an ordering or off-by-one mistake would
    live — under test without DESeq2 or a container.
    """

    @staticmethod
    def _extract(script_name: str) -> str:
        source = resources.files("txsuite.resources.bulk_r").joinpath(
            script_name
        ).read_text(encoding="utf-8")
        start = source.index("expanded_contrasts <- function")
        return source[start : source.index("\n}\n", start) + 3]

    def _pairs(self, script_name, mode, levels, reference, test) -> list[str]:
        harness = "\n".join(
            [
                self._extract(script_name),
                "pairs <- expanded_contrasts(%s, c(%s), %s, %s)"
                % (
                    _r_string(mode),
                    ", ".join(_r_string(level) for level in levels),
                    _r_string(reference),
                    _r_string(test),
                ),
                'cat(paste(sapply(pairs, function(p) paste0(p[[1]], "_vs_", p[[2]])),'
                ' collapse = "\\n"))',
            ]
        )
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "harness.R"
            script.write_text(harness, encoding="utf-8")
            completed = subprocess.run(
                ["Rscript", str(script)], capture_output=True, text=True, check=True
            )
        return [line for line in completed.stdout.splitlines() if line]

    def test_expansion_excludes_the_primary_contrast(self) -> None:
        levels = ("control", "low", "high")
        for script in ("deseq2.R", "alternative_de.R"):
            with self.subTest(script=script):
                self.assertEqual(
                    self._pairs(script, "single", levels, "control", "low"), []
                )
                self.assertEqual(
                    self._pairs(script, "vs-reference", levels, "control", "low"),
                    ["high_vs_control"],
                )
                self.assertEqual(
                    self._pairs(script, "all-pairs", levels, "control", "low"),
                    ["high_vs_control", "high_vs_low"],
                )

    def test_primary_choice_only_changes_which_pair_is_omitted(self) -> None:
        levels = ("control", "low", "high")
        self.assertEqual(
            self._pairs("deseq2.R", "all-pairs", levels, "control", "high"),
            ["low_vs_control", "high_vs_low"],
        )

    def test_expansion_is_correct_when_the_reference_sorts_last(self) -> None:
        """Expansion holds for any reference position, given ordered levels.

        This exercises the helper with levels ordered the way both callers now
        order them. It does not by itself prove the caller orders them: that is
        covered behaviourally by the mirrored-primary test below, which passes
        raw alphabetical levels, and textually by
        ``test_alternative_script_orders_levels_reference_first``.
        """

        levels = ("control", "treated", "untreated")
        for script in ("deseq2.R", "alternative_de.R"):
            for reference, test in (
                ("control", "treated"),
                ("treated", "control"),
                ("untreated", "control"),
            ):
                with self.subTest(script=script, reference=reference, test=test):
                    ordered = (reference,) + tuple(
                        level for level in levels if level != reference
                    )
                    pairs = self._pairs(script, "all-pairs", ordered, reference, test)
                    self.assertNotIn(f"{test}_vs_{reference}", pairs)
                    self.assertNotIn(f"{reference}_vs_{test}", pairs)
                    self.assertEqual(len(pairs), 2)
                    self.assertEqual(len(set(pairs)), 2)

    def test_mirrored_primary_is_dropped_even_from_unordered_levels(self) -> None:
        # Defence in depth: the filter itself rejects both orientations, so a
        # future caller that forgets to order the levels cannot reintroduce the
        # duplicate contrast.
        pairs = self._pairs(
            "alternative_de.R", "all-pairs", ("control", "treated", "untreated"),
            "treated", "control",
        )
        self.assertNotIn("treated_vs_control", pairs)
        self.assertNotIn("control_vs_treated", pairs)

    def test_two_level_designs_expand_to_nothing_beyond_the_primary(self) -> None:
        for mode in ("vs-reference", "all-pairs"):
            with self.subTest(mode=mode):
                self.assertEqual(
                    self._pairs("deseq2.R", mode, ("a", "b"), "a", "b"), []
                )

    def test_all_pairs_covers_every_combination(self) -> None:
        levels = tuple("abcde")
        pairs = self._pairs("deseq2.R", "all-pairs", levels, "a", "b")
        # Five levels give ten combinations; the primary is written separately.
        self.assertEqual(len(pairs), 9)
        self.assertEqual(len(set(pairs)), 9)
        self.assertNotIn("b_vs_a", pairs)


def _r_string(value: str) -> str:
    return '"%s"' % value


if __name__ == "__main__":
    unittest.main()
