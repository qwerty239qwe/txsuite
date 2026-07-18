from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def read_input(path: Path):
    if path.is_dir():
        if (path / "matrix.mtx").is_file():
            adata = sc.read_mtx(path / "matrix.mtx").T
            features = pd.read_csv(path / "features.tsv", sep="\t", header=None)
            barcodes = pd.read_csv(path / "barcodes.tsv", sep="\t", header=None)
            adata.var_names = features.iloc[:, 1].astype(str)
            adata.var["gene_ids"] = features.iloc[:, 0].astype(str).to_numpy()
            adata.obs_names = barcodes.iloc[:, 0].astype(str)
            return adata
        return sc.read_10x_mtx(path, var_names="gene_symbols")
    if path.suffix == ".h5ad":
        return sc.read_h5ad(path)
    if path.suffix == ".h5":
        return sc.read_10x_h5(path)
    raise ValueError("Input must be a 10x matrix directory, 10x .h5, or .h5ad")


def analyze(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    adata = read_input(Path(args.input))
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["pct_counts_mt"] <= args.max_mito_pct].copy()
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError("QC left fewer than 3 cells or 3 genes")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2_000, adata.n_vars))
    use_hvg = int(adata.var["highly_variable"].sum()) >= 3
    variables = int(adata.var["highly_variable"].sum()) if use_hvg else adata.n_vars
    n_comps = min(50, adata.n_obs - 1, variables - 1)
    if n_comps < 2:
        raise ValueError("Not enough cells or genes remain for PCA")
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=use_hvg)
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1), n_pcs=n_comps)
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    sc.tl.umap(adata, random_state=0)

    adata.write_h5ad(output / "analysis.h5ad", compression="gzip")
    adata.obs.to_csv(output / "cell-qc.tsv", sep="\t", index_label="barcode")
    clusters = pd.DataFrame(
        {
            "barcode": adata.obs_names,
            "leiden": adata.obs["leiden"].astype(str).to_numpy(),
            "umap_1": adata.obsm["X_umap"][:, 0],
            "umap_2": adata.obsm["X_umap"][:, 1],
        }
    )
    clusters.to_csv(output / "clusters.tsv", sep="\t", index=False)
    (output / "summary.json").write_text(
        json.dumps(
            {
                "scanpy": sc.__version__,
                "cells": adata.n_obs,
                "genes": adata.n_vars,
                "clusters": int(adata.obs["leiden"].nunique()),
            },
            indent=2,
        )
        + "\n"
    )


def pseudobulk(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    for column in (args.sample_column, args.design):
        if column not in adata.obs:
            raise ValueError(f"AnnData obs is missing column: {column}")
    if (
        adata.obs.groupby(args.sample_column, observed=True)[args.design].nunique() > 1
    ).any():
        raise ValueError("Each sample must have exactly one design value")
    matrix = adata.layers.get("counts", adata.X)
    samples = sorted(adata.obs[args.sample_column].astype(str).unique())
    columns = {}
    for sample in samples:
        mask = adata.obs[args.sample_column].astype(str).to_numpy() == sample
        columns[sample] = np.asarray(matrix[mask].sum(axis=0)).ravel()
    counts = pd.DataFrame(columns, index=adata.var_names)
    rounded = np.rint(counts.to_numpy())
    if (rounded < 0).any() or not np.allclose(counts.to_numpy(), rounded):
        raise ValueError("Pseudobulk requires non-negative integer raw counts")
    counts.iloc[:, :] = rounded.astype(np.int64)
    counts.to_csv(output / "pseudobulk-counts.tsv", sep="\t", index_label="gene_id")
    metadata = (
        adata.obs[[args.sample_column, args.design]]
        .assign(
            **{args.sample_column: lambda frame: frame[args.sample_column].astype(str)}
        )
        .drop_duplicates()
        .sort_values(args.sample_column)
        .rename(columns={args.sample_column: "sample"})
    )
    metadata.to_csv(output / "pseudobulk-metadata.tsv", sep="\t", index=False)


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(dest="command", required=True)
    analysis = commands.add_parser("analyze")
    analysis.add_argument("input")
    analysis.add_argument("output")
    analysis.add_argument("--min-genes", type=int, default=200)
    analysis.add_argument("--min-cells", type=int, default=3)
    analysis.add_argument("--max-mito-pct", type=float, default=20)
    analysis.add_argument("--resolution", type=float, default=1)
    bulk = commands.add_parser("pseudobulk")
    bulk.add_argument("input")
    bulk.add_argument("output")
    bulk.add_argument("--sample-column", required=True)
    bulk.add_argument("--design", required=True)
    return root


if __name__ == "__main__":
    args = parser().parse_args()
    {"analyze": analyze, "pseudobulk": pseudobulk}[args.command](args)
