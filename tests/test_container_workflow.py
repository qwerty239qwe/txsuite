from pathlib import Path
import unittest


WORKFLOW = Path(__file__).parents[1] / ".github" / "workflows" / "containers.yml"


class ContainerWorkflowTests(unittest.TestCase):
    def test_all_owned_images_are_published_with_digests(self) -> None:
        workflow = WORKFLOW.read_text(encoding="utf-8")

        for context in ("bulk_r", "single_cell_python", "spatial_python"):
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


if __name__ == "__main__":
    unittest.main()
