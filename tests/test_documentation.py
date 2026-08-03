from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class DocumentationTests(unittest.TestCase):
    def test_project_documentation_set_is_linked_from_readme(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        required = (
            "roadmap",
            "compatibility",
            "release",
            "agent-workflows",
            "workflow-schema",
            "results-layout",
        )
        for name in required:
            with self.subTest(document=name):
                path = ROOT / "docs" / f"{name}.md"
                self.assertTrue(path.is_file(), path)
                self.assertIn(f"docs/{name}.md", readme)

    def test_project_commands_and_resume_contract_are_documented(self) -> None:
        documents = "\n".join(
            (ROOT / relative).read_text(encoding="utf-8")
            for relative in ("README.md", "docs/agent-workflows.md")
        )
        commands = (
            "txsuite project init --preset",
            "txsuite project validate workflow.toml",
            "txsuite project plan workflow.toml",
            "txsuite project run workflow.toml",
            "txsuite project status",
            "--resume --run-id",
        )
        for command in commands:
            with self.subTest(command=command):
                self.assertIn(command, documents)

    def test_security_layout_and_release_blockers_are_explicit(self) -> None:
        layout = (ROOT / "docs" / "results-layout.md").read_text(encoding="utf-8")
        for name in (
            "workflow.json",
            "resolved-config.json",
            "command-plan.json",
            "run-manifest.json",
            "stage-results",
            "stdout.log",
            "stderr.log",
        ):
            self.assertIn(name, layout)
        self.assertIn("redact", layout.casefold())
        release = (ROOT / "docs" / "release.md").read_text(encoding="utf-8")
        for blocker in ("Docker", "Apptainer", "Space Ranger", "sha256"):
            self.assertIn(blocker, release)


if __name__ == "__main__":
    unittest.main()
