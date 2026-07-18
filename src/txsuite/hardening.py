from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from contextlib import closing
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit
from urllib.request import urlopen

from txsuite.runtime import TxSuiteError


SHA256 = re.compile(r"[0-9a-f]{64}")
LOCKED_IMAGE = re.compile(r"(?:@sha256:|^sha256:)[0-9a-f]{64}$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cache_reference(source: str, expected_sha256: str, name: str, root: Path) -> Path:
    expected = expected_sha256.lower()
    if not SHA256.fullmatch(expected):
        raise TxSuiteError(
            "Reference SHA-256 must be exactly 64 hexadecimal characters"
        )
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
        raise TxSuiteError(
            "Reference name may contain only letters, numbers, ., _ and -"
        )
    root.mkdir(parents=True, exist_ok=True)
    destination = root / name
    if destination.exists():
        actual = sha256_file(destination)
        if actual != expected:
            raise TxSuiteError(
                f"Cached reference checksum mismatch: expected {expected}, got {actual}"
            )
    else:
        with tempfile.TemporaryDirectory(dir=root) as directory:
            temporary = Path(directory) / name
            try:
                if "://" in source:
                    with (
                        closing(urlopen(source)) as response,
                        temporary.open("wb") as output,
                    ):
                        shutil.copyfileobj(response, output)
                else:
                    source_path = Path(source)
                    if not source_path.is_file():
                        raise TxSuiteError(f"Reference source does not exist: {source}")
                    shutil.copyfile(source_path, temporary)
            except (OSError, ValueError) as exc:
                raise TxSuiteError(f"Cannot cache reference: {exc}") from exc
            actual = sha256_file(temporary)
            if actual != expected:
                raise TxSuiteError(
                    f"Downloaded reference checksum mismatch: expected {expected}, got {actual}"
                )
            os.replace(temporary, destination)
    parsed = urlsplit(source)
    recorded_source = (
        urlunsplit((parsed.scheme, parsed.netloc, parsed.path, "", ""))
        if "://" in source
        else source
    )
    (root / f"{name}.json").write_text(
        json.dumps(
            {
                "name": name,
                "path": str(destination.resolve()),
                "sha256": expected,
                "source": recorded_source,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return destination


def image_is_locked(reference: str) -> bool:
    return bool(LOCKED_IMAGE.search(reference))
