from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

from txsuite.runtime import run_command
from txsuite.spatial import analysis_command, spacemake_command, spaceranger_command


class Phase3Test(unittest.TestCase):
    def test_spatial_commands_and_working_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            transcriptome = root / "reference"
            fastqs = root / "fastqs"
            project = root / "spacemake"
            for path in (transcriptome, fastqs, project):
                path.mkdir()
            image = root / "tissue.tif"
            image.touch()

            command = spaceranger_command(
                run_id="sample-1",
                transcriptome=transcriptome,
                fastqs=fastqs,
                image=image,
                sample="sample",
                slide="V19J01-123",
                area="A1",
                unknown_slide=False,
                create_bam=False,
                cores=8,
                memory=64,
            )
            self.assertEqual(command[:3], ["spaceranger", "count", "--id=sample-1"])
            self.assertIn("--create-bam=false", command)
            self.assertEqual(
                spacemake_command(project, 2),
                ["spacemake", "run", "--cores", "2", "--keep-going"],
            )

            zarr = root / "sample.zarr"
            zarr.mkdir()
            analyze = analysis_command(
                image="txsuite/spatial-python:test",
                input_path=zarr,
                outdir=root / "output",
                dataset_id="sample",
                min_counts=0,
                min_spots=0,
            )
            self.assertTrue(
                any("target=/input/data.zarr,readonly" in part for part in analyze)
            )

            run_dir = root / "run"
            run_command(
                [sys.executable, "-c", "from pathlib import Path; print(Path.cwd())"],
                run_dir=run_dir,
                task="cwd-test",
                backend="python",
                inputs={},
                outputs={},
                artifacts=[],
                cwd=project,
            )
            self.assertEqual(
                Path((run_dir / "stdout.log").read_text().strip()), project.resolve()
            )


if __name__ == "__main__":
    unittest.main()
