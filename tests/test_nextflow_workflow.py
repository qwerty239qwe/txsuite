from __future__ import annotations

import tempfile
import unittest
from importlib import resources
from pathlib import Path

from txsuite.config import DEFAULT_CONFIG
from txsuite.single_cell import pseudobulk_workflow_command


class NextflowWorkflowTests(unittest.TestCase):
    def test_packaged_dag_and_command(self) -> None:
        package = resources.files("txsuite.resources.nextflow")
        main = package.joinpath("main.nf").read_text(encoding="utf-8")
        pseudobulk = package.joinpath("modules/pseudobulk.nf").read_text(
            encoding="utf-8"
        )
        deseq2 = package.joinpath("modules/deseq2.nf").read_text(encoding="utf-8")
        config = package.joinpath("nextflow.config").read_text(encoding="utf-8")

        self.assertIn("PSEUDOBULK(input_ch)", main)
        self.assertIn("DESEQ2(PSEUDOBULK.out.counts", main)
        self.assertIn("container params.single_cell_image", pseudobulk)
        self.assertIn("container params.bulk_image", deseq2)
        self.assertIn("stub:", pseudobulk)
        self.assertIn("stub:", deseq2)
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
                nextflow_config=extra_config,
                resume=True,
            )

            self.assertEqual(command[:2], ["nextflow", "run"])
            self.assertIn("main.nf", command[2])
            self.assertIn("-resume", command)
            self.assertIn(str(extra_config.resolve()), command)
            self.assertIn(DEFAULT_CONFIG["images"]["single_cell_python"], command)
            self.assertIn(DEFAULT_CONFIG["images"]["bulk_r"], command)


if __name__ == "__main__":
    unittest.main()
