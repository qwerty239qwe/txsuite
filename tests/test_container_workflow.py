import unittest
from pathlib import Path

WORKFLOW = Path(__file__).parents[1] / ".github" / "workflows" / "containers.yml"
CI = Path(__file__).parents[1] / ".github" / "workflows" / "ci.yml"


class ContainerWorkflowTests(unittest.TestCase):
    def test_owned_images_accept_explicit_commands(self) -> None:
        root = Path(__file__).parents[1] / "src" / "txsuite" / "resources"
        for context in (
            "bulk_r",
            "biodbs",
            "salmon",
            "single_cell_python",
            "spatial_python",
            "star",
        ):
            dockerfile = (root / context / "Dockerfile").read_text(encoding="utf-8")
            self.assertNotIn("ENTRYPOINT", dockerfile)
            self.assertIn("WORKDIR /work", dockerfile)
            self.assertIn("CMD [", dockerfile)

    def test_all_owned_images_are_published_with_digests(self) -> None:
        workflow = WORKFLOW.read_text(encoding="utf-8")

        for context in (
            "bulk_r",
            "biodbs",
            "salmon",
            "single_cell_python",
            "spatial_python",
            "star",
        ):
            self.assertIn(f"src/txsuite/resources/{context}", workflow)
        for action in (
            "docker/setup-buildx-action@v4",
            "docker/login-action@v4",
            "docker/metadata-action@v6",
            "docker/build-push-action@v7",
        ):
            self.assertIn(action, workflow)
        self.assertIn("packages: write", workflow)
        self.assertIn("steps.build.outputs.digest", workflow)
        ci = CI.read_text(encoding="utf-8")
        for expected in (
            "single-cell-smoke:",
            "--doublets score",
            "--integration harmony",
            "X_pca_harmony",
            "marker-genes.tsv",
            "--group-column cell_type",
            "--covariate batch",
        ):
            self.assertIn(expected, ci)


if __name__ == "__main__":
    unittest.main()
