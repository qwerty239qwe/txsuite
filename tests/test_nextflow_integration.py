"""Run with TXSUITE_NEXTFLOW_TESTS=1 on a host with Nextflow and Python 3."""

import csv
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


@unittest.skipUnless(
    os.environ.get("TXSUITE_NEXTFLOW_TESTS") == "1", "Nextflow integration is opt-in"
)
class NextflowIntegrationTest(unittest.TestCase):
    def test_success_partial_failure_all_failed_and_resume(self):
        workflow = (
            Path(__file__).resolve().parents[1]
            / "src/txsuite/resources/nextflow/main.nf"
        )
        for failed in ((), ("b",), ("a", "b")):
            with (
                self.subTest(failed=failed),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                (root / "input.h5ad").touch()
                (root / "manifest.tsv").write_text(
                    "comparison\tdesign\treference\ttest\n"
                    "a\tcondition\tcontrol\ttreated\nb\tcondition\tcontrol\ttreated\n"
                )
                # Fault injection stays in test configuration, not production modules.
                names = ",".join(repr(name) for name in failed)
                (root / "test.config").write_text(
                    "process { withName: PSEUDOBULK { beforeScript = { "
                    f"[{names}].contains(task.tag) ? 'exit 1' : 'true'"
                    " } } }\n"
                )
                command = [
                    "nextflow",
                    "run",
                    str(workflow),
                    "-stub-run",
                    "-profile",
                    "local",
                    "-c",
                    str(root / "test.config"),
                    "-work-dir",
                    str(root / "work"),
                    "--input",
                    str(root / "input.h5ad"),
                    "--manifest",
                    str(root / "manifest.tsv"),
                    "--outdir",
                    str(root / "out"),
                    "--single_cell_image",
                    "stub",
                    "--bulk_image",
                    "stub",
                    "-with-trace",
                    str(root / "trace.tsv"),
                ]
                result = subprocess.run(
                    command,
                    cwd=root,
                    capture_output=True,
                    text=True,
                    timeout=180,
                    check=False,
                )
                self.assertEqual(
                    result.returncode == 0, not failed, result.stdout + result.stderr
                )
                with (root / "out/comparison-index.tsv").open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(
                    {r["comparison"]: r["status"] for r in rows},
                    {
                        name: "failed" if name in failed else "success"
                        for name in ("a", "b")
                    },
                )
                for row in rows:
                    if row["status"] == "success":
                        self.assertTrue((root / "out" / row["result"]).is_file())
                if not failed:
                    resumed = subprocess.run(
                        command[:-1] + [str(root / "resume-trace.tsv"), "-resume"],
                        cwd=root,
                        capture_output=True,
                        text=True,
                        timeout=180,
                        check=False,
                    )
                    self.assertEqual(
                        resumed.returncode, 0, resumed.stdout + resumed.stderr
                    )
                    with (root / "resume-trace.tsv").open() as handle:
                        trace = list(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual(len(trace), 5)
                    self.assertTrue(
                        all(row["status"] == "CACHED" for row in trace), trace
                    )
