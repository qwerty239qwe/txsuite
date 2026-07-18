from __future__ import annotations

import json
import shlex
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


class TxSuiteError(RuntimeError):
    pass


def format_command(command: list[str]) -> str:
    return shlex.join(command)


def run_command(
    command: list[str],
    *,
    run_dir: Path,
    task: str,
    backend: str,
    inputs: dict[str, Any],
    outputs: dict[str, Any],
    artifacts: list[dict[str, str]],
    cwd: Path | None = None,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    started_at = datetime.now(UTC).isoformat()
    base = {
        "schema_version": 1,
        "task": task,
        "backend": backend,
        "command": command,
        "inputs": inputs,
        "outputs": outputs,
        "started_at": started_at,
    }
    (run_dir / "command.txt").write_text(
        format_command(command) + "\n", encoding="utf-8"
    )
    _write_json(run_dir / "run.json", {**base, "status": "running"})

    try:
        with (
            (run_dir / "stdout.log").open("w", encoding="utf-8") as stdout,
            (run_dir / "stderr.log").open("w", encoding="utf-8") as stderr,
        ):
            result = subprocess.run(
                command,
                stdout=stdout,
                stderr=stderr,
                text=True,
                check=False,
                cwd=cwd,
            )
    except OSError as exc:
        _write_json(
            run_dir / "run.json",
            {
                **base,
                "status": "error",
                "duration_seconds": time.time() - started,
                "error": str(exc),
            },
        )
        raise TxSuiteError(f"Cannot run {command[0]}: {exc}") from exc

    finished = {
        **base,
        "status": "success" if result.returncode == 0 else "failed",
        "exit_code": result.returncode,
        "duration_seconds": time.time() - started,
    }
    _write_json(run_dir / "run.json", finished)
    if result.returncode:
        raise TxSuiteError(
            f"{task} failed with exit code {result.returncode}; see {run_dir / 'stderr.log'}"
        )
    _write_json(
        run_dir / "txsuite-results.json",
        {
            "schema_version": 1,
            "task": task,
            "backend": backend,
            "artifacts": artifacts,
            "run": "run.json",
        },
    )


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
