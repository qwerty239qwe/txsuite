from __future__ import annotations

import argparse
import json
from pathlib import Path

import scanpy as sc
import spatialdata
import spatialdata_io
import squidpy as sq
from spatialdata_io.experimental import to_legacy_anndata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--dataset-id", default="sample")
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--min-spots", type=int, default=3)
    args = parser.parse_args()

    source = Path(args.input)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    if source.name.endswith(".zarr"):
        sdata = spatialdata.read_zarr(source)
    else:
        sdata = spatialdata_io.visium(source, dataset_id=args.dataset_id)
    sdata.write(output / "spatialdata.zarr")

    coordinate_system = (
        "global"
        if "global" in sdata.coordinate_systems
        else next(iter(sdata.coordinate_systems))
    )
    adata = to_legacy_anndata(sdata, coordinate_system=coordinate_system)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_genes(adata, min_cells=args.min_spots)
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError("QC left fewer than 3 spots or 3 genes")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10_000)
    sc.pp.log1p(adata)
    sq.gr.spatial_neighbors_grid(adata, n_neighs=min(6, adata.n_obs - 1))

    adata.write_h5ad(output / "analysis.h5ad", compression="gzip")
    adata.obs.to_csv(output / "spot-qc.tsv", sep="\t", index_label="barcode")
    (output / "summary.json").write_text(
        json.dumps(
            {
                "spatialdata": spatialdata.__version__,
                "spatialdata_io": spatialdata_io.__version__,
                "squidpy": sq.__version__,
                "spots": adata.n_obs,
                "genes": adata.n_vars,
                "spatial_edges": int(adata.obsp["spatial_connectivities"].nnz),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
