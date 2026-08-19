"""Convert an alevin-fry quantification directory into an AnnData ``.h5ad``.

``simpleaf quant`` writes a matrix-market bundle rather than the ``.h5ad`` that
every downstream TxSuite single-cell stage consumes. The format is small and
stable, so it is parsed directly here instead of adding another dependency to
the image.

USA-mode output stores spliced, unspliced, and ambiguous counts as three column
blocks over the same gene order. ``X`` is built as spliced + ambiguous, the
convention alevin-fry documents for standard single-cell gene expression, and
all three blocks are kept as layers so a caller can re-derive any other choice.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.io
import scipy.sparse as sp


class ConversionError(RuntimeError):
    """The alevin-fry output directory is missing files or internally inconsistent."""


def _quant_root(path: Path) -> Path:
    """Accept a simpleaf run directory, its af_quant dir, or the alevin dir itself."""

    for candidate in (path / "af_quant" / "alevin", path / "alevin", path):
        if (candidate / "quants_mat.mtx").is_file():
            return candidate
    raise ConversionError(f"No quants_mat.mtx found under {path}")


def _read_lines(path: Path) -> list[str]:
    if not path.is_file():
        raise ConversionError(f"Missing alevin-fry file: {path}")
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def load_alevin(path: Path) -> ad.AnnData:
    """Load an alevin-fry quantification into an AnnData object."""

    root = _quant_root(path)
    barcodes = _read_lines(root / "quants_mat_rows.txt")
    columns = _read_lines(root / "quants_mat_cols.txt")
    matrix = sp.csr_matrix(scipy.io.mmread(root / "quants_mat.mtx"))

    if matrix.shape != (len(barcodes), len(columns)):
        raise ConversionError(
            f"Matrix shape {matrix.shape} does not match "
            f"{len(barcodes)} barcodes by {len(columns)} columns"
        )

    metadata_path = root / "quants_mat.json"
    usa_mode = False
    if metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        usa_mode = bool(metadata.get("usa_mode", False))

    layers: dict[str, sp.spmatrix] = {}
    if usa_mode:
        if len(columns) % 3 != 0:
            raise ConversionError(
                "USA-mode output must have a column count divisible by three, "
                f"got {len(columns)}"
            )
        width = len(columns) // 3
        spliced = matrix[:, :width]
        unspliced = matrix[:, width : 2 * width]
        ambiguous = matrix[:, 2 * width :]
        genes = [name.rsplit("-", 1)[0] if name.endswith("-S") else name for name in columns[:width]]
        counts = (spliced + ambiguous).tocsr()
        layers = {
            "spliced": spliced.tocsr(),
            "unspliced": unspliced.tocsr(),
            "ambiguous": ambiguous.tocsr(),
        }
    else:
        genes = columns
        counts = matrix

    counts = counts.astype(np.float32)
    adata = ad.AnnData(X=counts)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    adata.layers["counts"] = counts.copy()
    for name, layer in layers.items():
        adata.layers[name] = layer.astype(np.float32)
    adata.uns["alevin"] = {"usa_mode": usa_mode, "source": str(root)}
    return adata


def load_samples(quant_dirs: list[Path]) -> ad.AnnData:
    """Load one or more quantifications into a single sample-labelled AnnData.

    The sample name is the quantification directory name, which is how both the
    native DAG and ``simpleaf quant --output`` lay results out. Barcodes are
    prefixed with the sample so a concatenated object keeps unique observations.
    """

    if not quant_dirs:
        raise ConversionError("At least one quantification directory is required")

    loaded: dict[str, ad.AnnData] = {}
    for quant_dir in quant_dirs:
        sample = quant_dir.name
        if sample in loaded:
            raise ConversionError(f"Duplicate sample name: {sample}")
        adata = load_alevin(quant_dir)
        adata.obs["sample"] = sample
        adata.obs["sample"] = adata.obs["sample"].astype("category")
        loaded[sample] = adata

    if len(loaded) == 1:
        return next(iter(loaded.values()))
    return ad.concat(
        loaded,
        axis=0,
        join="outer",
        label=None,
        index_unique="-",
        fill_value=0,
        merge="first",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "quant_dirs",
        nargs="+",
        type=Path,
        help="simpleaf/alevin-fry output directories; the directory name is the sample",
    )
    parser.add_argument("--output", required=True, type=Path, help="destination .h5ad path")
    arguments = parser.parse_args(argv)
    try:
        adata = load_samples(arguments.quant_dirs)
    except ConversionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
