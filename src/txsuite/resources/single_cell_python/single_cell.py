from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


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


def add_metadata(adata, path: Path, barcode_column: str) -> list[str]:
    metadata = pd.read_csv(path, sep="\t")
    if barcode_column not in metadata:
        raise ValueError(f"Metadata is missing barcode column: {barcode_column}")
    barcodes = metadata[barcode_column]
    if barcodes.isna().any() or barcodes.astype(str).str.strip().eq("").any():
        raise ValueError("Metadata barcodes must be non-empty")
    barcodes = barcodes.astype(str)
    if barcodes.duplicated().any():
        raise ValueError("Metadata barcodes must be unique")
    cells = pd.Index(adata.obs_names.astype(str))
    metadata_cells = pd.Index(barcodes)
    missing = cells.difference(metadata_cells)
    extra = metadata_cells.difference(cells)
    if len(missing) or len(extra):
        raise ValueError(
            "Metadata barcodes must exactly match AnnData cells "
            f"(missing={len(missing)}, extra={len(extra)})"
        )
    columns = [column for column in metadata if column != barcode_column]
    overlap = sorted(set(columns) & set(adata.obs.columns))
    if overlap:
        raise ValueError(
            f"Metadata would replace AnnData columns: {', '.join(overlap)}"
        )
    metadata = metadata.assign(**{barcode_column: barcodes}).set_index(barcode_column)
    metadata = metadata.loc[cells]
    for column in columns:
        adata.obs[column] = metadata[column].to_numpy()
    return columns


def prepare_counts(adata, counts_layer: str) -> str:
    if counts_layer not in adata.layers and counts_layer != "counts":
        raise ValueError(f"Counts layer does not exist: {counts_layer}")
    source = f"layer:{counts_layer}" if counts_layer in adata.layers else "X"
    counts = (
        adata.layers[counts_layer].copy()
        if counts_layer in adata.layers
        else adata.X.copy()
    )
    validate_counts(counts, source)
    adata.X = counts.copy()
    adata.layers["counts"] = counts
    return source


def validate_counts(counts, source: str) -> None:
    values = counts.data if sparse.issparse(counts) else np.asarray(counts)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.equal(values, np.rint(values)).all()
        or (values >= 2**63).any()
    ):
        raise ValueError(
            f"{source} must contain non-negative integer raw counts; "
            "provide --counts-layer for normalized H5AD input"
        )


def reset_analysis(adata) -> None:
    for key in ("leiden", "doublet_score", "predicted_doublet"):
        if key in adata.obs:
            del adata.obs[key]
    for key in (
        "highly_variable",
        "highly_variable_intersection",
        "highly_variable_nbatches",
        "means",
        "dispersions",
        "dispersions_norm",
    ):
        if key in adata.var:
            del adata.var[key]
    for key in ("X_pca", "X_pca_harmony", "X_umap"):
        adata.obsm.pop(key, None)
    adata.varm.pop("PCs", None)
    for key in ("distances", "connectivities"):
        adata.obsp.pop(key, None)
    for key in (
        "hvg",
        "leiden",
        "log1p",
        "neighbors",
        "pca",
        "rank_genes_groups",
        "scrublet",
        "umap",
    ):
        adata.uns.pop(key, None)


def analyze(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    if not 0 < args.expected_doublet_rate < 1:
        raise ValueError("Expected doublet rate must be in (0, 1)")
    if args.doublet_threshold is not None and (
        not np.isfinite(args.doublet_threshold) or args.doublet_threshold < 0
    ):
        raise ValueError("Doublet threshold must be non-negative")
    if args.top_markers < 1:
        raise ValueError("Top markers must be positive")
    if args.integration == "harmony" and not args.batch_column:
        raise ValueError("Harmony integration requires a batch column")
    if args.integration == "harmony" and args.stop_after == "qc":
        raise ValueError("Harmony integration requires PCA; stop-after cannot be qc")
    if args.target_sum <= 0 or args.n_hvg < 3 or args.n_pcs < 2 or args.n_neighbors < 2:
        raise ValueError(
            "Target sum must be positive; HVGs >= 3; PCs and neighbors >= 2"
        )
    if not 0 <= args.umap_min_dist <= 1:
        raise ValueError("UMAP minimum distance must be in [0, 1]")
    adata = read_input(Path(args.input))
    adata.var_names_make_unique()
    counts_source = prepare_counts(adata, args.counts_layer)
    reset_analysis(adata)
    metadata_columns = (
        add_metadata(adata, Path(args.metadata), args.barcode_column)
        if args.metadata
        else []
    )
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["pct_counts_mt"] <= args.max_mito_pct].copy()
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError("QC left fewer than 3 cells or 3 genes")

    doublets = 0
    if args.doublets != "off":
        if args.doublet_batch_column and args.doublet_batch_column not in adata.obs:
            raise ValueError(
                "AnnData obs is missing doublet batch column: "
                f"{args.doublet_batch_column}"
            )
        if args.doublet_batch_column:
            doublet_batches = adata.obs[args.doublet_batch_column]
            if (
                doublet_batches.isna().any()
                or doublet_batches.astype(str).str.strip().eq("").any()
            ):
                raise ValueError("Doublet batch values must be non-empty")
        batch_limit = (
            int(doublet_batches.value_counts().min()) - 2
            if args.doublet_batch_column
            else adata.n_obs - 2
        )
        n_prin_comps = min(30, batch_limit, adata.n_vars - 2)
        if n_prin_comps < 1:
            raise ValueError("Scrublet requires at least 3 cells and 3 genes per batch")
        scrublet_args = {
            "batch_key": args.doublet_batch_column,
            "expected_doublet_rate": args.expected_doublet_rate,
            "n_prin_comps": n_prin_comps,
            "random_state": 0,
        }
        if args.doublet_threshold is not None:
            scrublet_args["threshold"] = args.doublet_threshold
        sc.pp.scrublet(adata, **scrublet_args)
        doublets = int(adata.obs["predicted_doublet"].sum())
        if args.doublets == "filter":
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            if adata.n_obs < 3:
                raise ValueError("Doublet filtering left fewer than 3 cells")

    batch_count = 0
    if args.batch_column:
        if args.batch_column not in adata.obs:
            raise ValueError(
                f"AnnData obs is missing batch column: {args.batch_column}"
            )
        batch_values = adata.obs[args.batch_column]
        if (
            batch_values.isna().any()
            or batch_values.astype(str).str.strip().eq("").any()
        ):
            raise ValueError("Batch values must be non-empty")
        batch_count = int(batch_values.nunique())
        if batch_count < 2:
            raise ValueError("Batch column must contain at least two batches")

    marker_columns = [
        "cluster",
        "gene",
        "score",
        "log2fc",
        "pvalue",
        "padj",
        "pct_in",
        "pct_out",
    ]
    markers = pd.DataFrame(columns=marker_columns)
    integration_basis = None
    has_clusters = False
    has_umap = False

    if args.stop_after != "qc":
        sc.pp.normalize_total(adata, target_sum=args.target_sum)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=min(args.n_hvg, adata.n_vars),
            flavor=args.hvg_flavor,
            batch_key=args.batch_column,
        )
        use_hvg = int(adata.var["highly_variable"].sum()) >= 3
        variables = int(adata.var["highly_variable"].sum()) if use_hvg else adata.n_vars
        n_comps = min(args.n_pcs, adata.n_obs - 1, variables - 1)
        if n_comps < 2:
            raise ValueError("Not enough cells or genes remain for PCA")
        sc.tl.pca(
            adata,
            n_comps=n_comps,
            mask_var="highly_variable" if use_hvg else None,
        )
        integration_basis = "X_pca"
        if args.integration == "harmony":
            import harmonypy

            # ponytail: Scanpy 1.12.2 transposes harmonypy>=0.1 output; remove
            # this direct call when pinned Scanpy has native Harmony support.
            result = harmonypy.run_harmony(
                adata.obsm["X_pca"],
                adata.obs,
                args.batch_column,
                nclust=max(2, min(round(adata.n_obs / 30), 100)),
                random_state=0,
                verbose=False,
            )
            corrected = np.asarray(result.Z_corr)
            if (
                corrected.shape != adata.obsm["X_pca"].shape
                or not np.isfinite(corrected).all()
            ):
                raise ValueError("Harmony returned invalid corrected PCs")
            integration_basis = "X_pca_harmony"
            adata.obsm[integration_basis] = corrected

        if args.stop_after in {"clusters", "all"}:
            neighbor_args = {"n_neighbors": min(args.n_neighbors, adata.n_obs - 1)}
            if integration_basis == "X_pca_harmony":
                neighbor_args["use_rep"] = integration_basis
            else:
                neighbor_args["n_pcs"] = n_comps
            sc.pp.neighbors(adata, **neighbor_args)
            sc.tl.leiden(
                adata,
                resolution=args.resolution,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            has_clusters = True

            if args.stop_after == "all" and not args.skip_umap:
                sc.tl.umap(adata, min_dist=args.umap_min_dist, random_state=0)
                has_umap = True

            if args.stop_after == "all" and not args.skip_markers:
                cluster_sizes = adata.obs["leiden"].value_counts()
                marker_groups = cluster_sizes[cluster_sizes >= 2].index.tolist()
                if adata.obs["leiden"].nunique() > 1 and marker_groups:
                    sc.tl.rank_genes_groups(
                        adata,
                        groupby="leiden",
                        groups=marker_groups,
                        method=args.marker_method,
                        n_genes=min(args.top_markers, adata.n_vars),
                        pts=True,
                        use_raw=False,
                    )
                    markers = sc.get.rank_genes_groups_df(adata, group=None).rename(
                        columns={
                            "group": "cluster",
                            "names": "gene",
                            "scores": "score",
                            "logfoldchanges": "log2fc",
                            "pvals": "pvalue",
                            "pvals_adj": "padj",
                            "pct_nz_group": "pct_in",
                            "pct_nz_reference": "pct_out",
                        }
                    )
                    markers = markers.reindex(columns=marker_columns)
    markers.to_csv(output / "marker-genes.tsv", sep="\t", index=False)

    adata.write_h5ad(output / "analysis.h5ad", compression="gzip")
    adata.obs.to_csv(output / "cell-qc.tsv", sep="\t", index_label="barcode")
    clusters = pd.DataFrame(columns=["barcode", "leiden", "umap_1", "umap_2"])
    if has_clusters:
        clusters = pd.DataFrame(
            {
                "barcode": adata.obs_names,
                "leiden": adata.obs["leiden"].astype(str).to_numpy(),
                "umap_1": adata.obsm["X_umap"][:, 0] if has_umap else np.nan,
                "umap_2": adata.obsm["X_umap"][:, 1] if has_umap else np.nan,
            }
        )
    clusters.to_csv(output / "clusters.tsv", sep="\t", index=False)
    (output / "summary.json").write_text(
        json.dumps(
            {
                "scanpy": version("scanpy"),
                "cells": adata.n_obs,
                "genes": adata.n_vars,
                "clusters": int(adata.obs["leiden"].nunique())
                if has_clusters
                else None,
                "metadata_columns": metadata_columns,
                "counts_source": counts_source,
                "batch_column": args.batch_column,
                "batches": batch_count,
                "integration": args.integration,
                "integration_basis": integration_basis,
                "harmonypy": (
                    version("harmonypy") if args.integration == "harmony" else None
                ),
                "doublet_mode": args.doublets,
                "predicted_doublets": doublets,
                "markers": len(markers),
                "marker_method": args.marker_method,
                "stop_after": args.stop_after,
                "umap": has_umap,
                "target_sum": args.target_sum,
                "n_hvg": args.n_hvg,
                "hvg_flavor": args.hvg_flavor,
                "n_pcs": args.n_pcs,
                "n_neighbors": args.n_neighbors,
                "umap_min_dist": args.umap_min_dist,
            },
            indent=2,
        )
        + "\n"
    )


def pseudobulk(args: argparse.Namespace) -> None:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(args.input)
    if bool(args.group_column) != bool(args.group_value):
        raise ValueError("Group column and group value must be used together")
    metadata_columns = [args.sample_column, args.design, *args.covariate]
    if args.group_column:
        metadata_columns.append(args.group_column)
    if len(set(metadata_columns)) != len(metadata_columns):
        raise ValueError("Sample, design, group, and covariate columns must differ")
    if bool(args.reference) != bool(args.test) or (
        args.reference is not None and args.reference == args.test
    ):
        raise ValueError("Reference and test must be used together and differ")
    for column in metadata_columns:
        if column not in adata.obs:
            raise ValueError(f"AnnData obs is missing column: {column}")
        values = adata.obs[column]
        if values.isna().any() or values.astype(str).str.strip().eq("").any():
            raise ValueError(f"AnnData obs column contains empty values: {column}")
    if args.group_column:
        adata = adata[
            adata.obs[args.group_column].astype(str) == args.group_value
        ].copy()
        if not adata.n_obs:
            raise ValueError(f"No cells match {args.group_column}={args.group_value}")
    fixed_columns = [args.design, *args.covariate]
    if args.group_column:
        fixed_columns.append(args.group_column)
    for column in fixed_columns:
        if (
            adata.obs.groupby(args.sample_column, observed=True)[column].nunique() > 1
        ).any():
            raise ValueError(f"Each sample must have exactly one {column} value")
    if args.counts_layer not in adata.layers and args.counts_layer != "counts":
        raise ValueError(f"Counts layer does not exist: {args.counts_layer}")
    matrix = adata.layers.get(args.counts_layer, adata.X)
    validate_counts(matrix, args.counts_layer)
    samples = sorted(adata.obs[args.sample_column].astype(str).unique())
    count_columns = {}
    for sample in samples:
        mask = adata.obs[args.sample_column].astype(str).to_numpy() == sample
        count_columns[sample] = np.asarray(matrix[mask].astype(np.float64).sum(axis=0)).ravel()
    counts = pd.DataFrame(count_columns, index=adata.var_names)
    # Float64 aggregation is exact only below this integer limit.
    if (counts.to_numpy() >= 2**53).any():
        raise ValueError("Pseudobulk totals exceed the exact integer range")
    validate_counts(counts.to_numpy(), "Pseudobulk totals")
    counts = counts.astype(np.int64)
    counts.to_csv(output / "pseudobulk-counts.tsv", sep="\t", index_label="gene_id")
    metadata = (
        adata.obs[metadata_columns]
        .assign(
            **{args.sample_column: lambda frame: frame[args.sample_column].astype(str)}
        )
        .drop_duplicates()
        .sort_values(args.sample_column)
        .rename(columns={args.sample_column: "sample"})
    )
    if args.reference and args.test:
        levels = metadata[args.design].astype(str).value_counts()
        for level in (args.reference, args.test):
            if levels.get(level, 0) < 2:
                raise ValueError(
                    f"Design level {level} requires at least two biological samples"
                )
    metadata.to_csv(output / "pseudobulk-metadata.tsv", sep="\t", index=False)


def collect_de(args: argparse.Namespace) -> None:
    from collect_de import collect_de as collect

    collect(args)


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
    analysis.add_argument("--metadata")
    analysis.add_argument("--barcode-column", default="barcode")
    analysis.add_argument("--batch-column")
    analysis.add_argument("--integration", choices=("none", "harmony"), default="none")
    analysis.add_argument("--counts-layer", default="counts")
    analysis.add_argument("--target-sum", type=float, default=10_000)
    analysis.add_argument("--n-hvg", type=int, default=2_000)
    analysis.add_argument(
        "--hvg-flavor", choices=("seurat", "cell_ranger"), default="seurat"
    )
    analysis.add_argument("--n-pcs", type=int, default=50)
    analysis.add_argument("--n-neighbors", type=int, default=15)
    analysis.add_argument("--umap-min-dist", type=float, default=0.5)
    analysis.add_argument(
        "--marker-method", choices=("wilcoxon", "t-test"), default="wilcoxon"
    )
    analysis.add_argument(
        "--stop-after", choices=("qc", "pca", "clusters", "all"), default="all"
    )
    analysis.add_argument("--skip-umap", action="store_true")
    analysis.add_argument("--skip-markers", action="store_true")
    analysis.add_argument(
        "--doublets", choices=("off", "score", "filter"), default="off"
    )
    analysis.add_argument("--doublet-batch-column")
    analysis.add_argument("--expected-doublet-rate", type=float, default=0.05)
    analysis.add_argument("--doublet-threshold", type=float)
    analysis.add_argument("--top-markers", type=int, default=100)
    bulk = commands.add_parser("pseudobulk")
    bulk.add_argument("input")
    bulk.add_argument("output")
    bulk.add_argument("--sample-column", required=True)
    bulk.add_argument("--design", required=True)
    bulk.add_argument("--counts-layer", default="counts")
    bulk.add_argument("--group-column")
    bulk.add_argument("--group-value")
    bulk.add_argument("--covariate", action="append", default=[])
    bulk.add_argument("--reference")
    bulk.add_argument("--test")
    collect = commands.add_parser("collect-de")
    collect.add_argument("expected")
    collect.add_argument("output")
    collect.add_argument("result_dir", nargs="*")
    return root


if __name__ == "__main__":
    args = parser().parse_args()
    {"analyze": analyze, "pseudobulk": pseudobulk, "collect-de": collect_de}[
        args.command
    ](args)
