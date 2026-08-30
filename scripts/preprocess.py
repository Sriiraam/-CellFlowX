from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]

INPUT = (
    ROOT
    / "data"
    / "processed"
    / "cellflowx_qc.h5ad"
)

OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "cellflowx_preprocessed.h5ad"
)


def main():

    print(f"Loading: {INPUT}")

    adata = sc.read_h5ad(INPUT)

    if adata.X is None:
        raise ValueError("adata.X is missing")

    if not sparse.issparse(adata.X):
        raise ValueError("Expected sparse count matrix")

    print(
        f"Input: {adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    # --------------------------------------------------
    # Preserve original raw counts
    # --------------------------------------------------

    adata.layers["counts"] = adata.X.copy()

    # --------------------------------------------------
    # Gene filtering
    # Require expression in at least 3 cells
    # --------------------------------------------------

    sc.pp.filter_genes(
        adata,
        min_cells=3,
    )

    print(
        f"After gene filtering: "
        f"{adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    # --------------------------------------------------
    # Library-size normalization
    # --------------------------------------------------

    sc.pp.normalize_total(
        adata,
        target_sum=1e4,
    )

    # --------------------------------------------------
    # Log transformation
    # --------------------------------------------------

    sc.pp.log1p(adata)

    # --------------------------------------------------
    # Highly variable genes
    #
    # Seurat flavor works on log-normalized data
    # and avoids needing scikit-misc.
    # --------------------------------------------------

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat",
        batch_key="geo_accession",
        subset=False,
    )

    n_hvg = int(
        adata.var["highly_variable"].sum()
    )

    print(
        f"Highly variable genes: {n_hvg:,}"
    )

    # --------------------------------------------------
    # Validation
    # --------------------------------------------------

    if "counts" not in adata.layers:
        raise ValueError("Raw counts layer missing")

    counts = adata.layers["counts"]

    if not sparse.issparse(counts):
        raise ValueError(
            "Counts layer is no longer sparse"
        )

    if not np.allclose(
        counts.data,
        np.round(counts.data),
    ):
        raise ValueError(
            "Counts layer is not integer-valued"
        )

    if "highly_variable" not in adata.var.columns:
        raise ValueError(
            "HVG annotation missing"
        )

    if n_hvg == 0:
        raise ValueError(
            "No HVGs were identified"
        )

    # --------------------------------------------------
    # Save preprocessing metadata
    # --------------------------------------------------

    adata.uns["preprocessing"] = {
        "gene_filter": "min_cells=3",
        "normalization": "normalize_total target_sum=10000",
        "transformation": "log1p",
        "hvg_method": "seurat",
        "hvg_n_top_genes": 2000,
        "hvg_batch_key": "geo_accession",
        "raw_counts_preserved": "layers['counts']",
    }

    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    adata.write_h5ad(
        OUTPUT,
        compression="gzip",
    )

    print("\n=== PREPROCESSING COMPLETE ===")
    print(
        f"Cells:              {adata.n_obs:,}"
    )
    print(
        f"Genes:              {adata.n_vars:,}"
    )
    print(
        f"HVGs:               {n_hvg:,}"
    )
    print(
        f"Counts layer sparse: "
        f"{sparse.issparse(adata.layers['counts'])}"
    )
    print(
        "Raw counts integer: True"
    )
    print(f"\nSaved: {OUTPUT}")


if __name__ == "__main__":
    main()
