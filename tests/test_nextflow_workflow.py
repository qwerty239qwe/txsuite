from __future__ import annotations

import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG
from txsuite.single_cell import (
    pseudobulk_manifest_workflow_command,
    pseudobulk_workflow_command,
)


class NextflowWorkflowTests(unittest.TestCase):
    def test_packaged_dag_and_command(self) -> None:
        package = resources.files("txsuite.resources.nextflow")
        main = package.joinpath("main.nf").read_text(encoding="utf-8")
        pseudobulk = package.joinpath("modules/pseudobulk.nf").read_text(
            encoding="utf-8"
        )
        bulk_de = package.joinpath("modules/bulk_de.nf").read_text(encoding="utf-8")
        collect = package.joinpath("modules/collect_de.nf").read_text(encoding="utf-8")
        subworkflow = package.joinpath("subworkflows/pseudobulk_de.nf").read_text(
            encoding="utf-8"
        )
        config = package.joinpath("nextflow.config").read_text(encoding="utf-8")

        self.assertIn("PSEUDOBULK_DE(comparisons)", main)
        self.assertIn("COLLECT_DE(expected, result_inputs, file(", main)
        self.assertIn("PSEUDOBULK(comparisons)", subworkflow)
        self.assertIn("BULK_DE(PSEUDOBULK.out.data)", subworkflow)
        self.assertNotIn("params.", pseudobulk)
        self.assertNotIn("params.", bulk_de)
        self.assertNotIn("params.", collect)
        self.assertNotIn("publishDir", pseudobulk)
        self.assertIn("--group-column", pseudobulk)
        self.assertIn("--covariate", pseudobulk)
        self.assertIn("alternative_de.R", bulk_de)
        self.assertIn("meta.method", bulk_de)
        self.assertIn('${collector}', collect)
        self.assertIn("stub:", pseudobulk)
        self.assertIn("stub:", bulk_de)
        self.assertIn("workflow.failOnIgnore = true", config)
        for profile in ("local", "docker", "apptainer"):
            self.assertIn(f"{profile} {{", config)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "data.h5ad"
            h5ad.touch()
            extra_config = root / "slurm.config"
            extra_config.write_text("process.executor = 'slurm'\n", encoding="utf-8")
            command = pseudobulk_workflow_command(
                DEFAULT_CONFIG,
                h5ad=h5ad,
                outdir=root / "results",
                sample_column="sample",
                design="condition",
                reference="control",
                test="treated",
                group_column="cell_type",
                group_value="T_cell",
                covariates=("batch",),
                nextflow_config=extra_config,
                resume=True,
            )

            self.assertEqual(command[:2], ["nextflow", "run"])
            self.assertIn("main.nf", command[2])
            self.assertIn("-resume", command)
            self.assertIn(str(extra_config.resolve()), command)
            self.assertIn(DEFAULT_CONFIG["images"]["single_cell_python"], command)
            self.assertIn(DEFAULT_CONFIG["images"]["bulk_r"], command)
            self.assertIn("cell_type", command)
            self.assertIn("T_cell", command)
            self.assertIn("batch", command)

            manifest = root / "comparisons.tsv"
            manifest.write_text(
                "comparison\tdesign\treference\ttest\tmethod\n"
                "t_cells\tcondition\tcontrol\ttreated\tedger\n",
                encoding="utf-8",
            )
            batch = pseudobulk_manifest_workflow_command(
                DEFAULT_CONFIG,
                manifest=manifest,
                h5ad=h5ad,
                outdir=root / "batch",
                sample_column="sample",
            )
            self.assertEqual(batch[:2], ["nextflow", "run"])
            self.assertIn("--manifest", batch)
            self.assertIn(str(manifest.resolve()), batch)


if __name__ == "__main__":
    unittest.main()
