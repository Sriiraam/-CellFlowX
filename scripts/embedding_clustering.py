import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    parser.add_argument("--pca-variance", required=True)
    parser.add_argument("--cluster-summary", required=True)
    parser.add_argument("--sample-cluster-summary", required=True)

    parser.add_argument("--pca-plot", required=True)
    parser.add_argument("--umap-sample", required=True)
    parser.add_argument("--umap-cluster", required=True)

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=15,
    )

    parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
    )

    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading: {args.input}")

    adata = sc.read_h5ad(args.input)

    print(
        f"Input: {adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    # --------------------------------------------------
    # VALIDATION
    # --------------------------------------------------

    if adata.X is None:
        raise ValueError("adata.X is missing")

    if not sparse.issparse(adata.X):
        raise ValueError(
            "Expected sparse normalized expression matrix"
        )

    if "highly_variable" not in adata.var.columns:
        raise ValueError(
            "highly_variable annotation missing"
        )

    n_hvg = int(
        adata.var["highly_variable"].sum()
    )

    print(f"HVGs available: {n_hvg:,}")

    if n_hvg < args.n_pcs:
        raise ValueError(
            "Number of HVGs is smaller than requested PCs"
        )

    if "counts" not in adata.layers:
        raise ValueError(
            "Raw counts layer missing"
        )

    if not sparse.issparse(
        adata.layers["counts"]
    ):
        raise ValueError(
            "Raw counts layer is not sparse"
        )

    # --------------------------------------------------
    # PCA
    # --------------------------------------------------

    print("\nRunning PCA...")

    sc.tl.pca(
        adata,
        n_comps=50,
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=42,
    )

    if "X_pca" not in adata.obsm:
        raise ValueError(
            "PCA embedding was not generated"
        )

    print(
        "PCA shape:",
        adata.obsm["X_pca"].shape,
    )

    # --------------------------------------------------
    # PCA VARIANCE
    # --------------------------------------------------

    variance_ratio = (
        adata.uns["pca"]["variance_ratio"]
    )

    pca_df = pd.DataFrame({
        "PC": np.arange(
            1,
            len(variance_ratio) + 1
        ),
        "variance_ratio": variance_ratio,
        "cumulative_variance":
            np.cumsum(variance_ratio),
    })

    Path(args.pca_variance).parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pca_df.to_csv(
        args.pca_variance,
        index=False,
    )

    # --------------------------------------------------
    # PCA VARIANCE PLOT
    # --------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.plot(
        pca_df["PC"],
        pca_df["cumulative_variance"],
        marker="o",
        markersize=3,
    )

    ax.set_xlabel(
        "Principal Component"
    )

    ax.set_ylabel(
        "Cumulative explained variance"
    )

    ax.set_title(
        "CellFlowX PCA Explained Variance"
    )

    fig.tight_layout()

    fig.savefig(
        args.pca_plot,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    # --------------------------------------------------
    # NEIGHBOR GRAPH
    # --------------------------------------------------

    print("\nBuilding neighbor graph...")

    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        random_state=42,
    )

    # --------------------------------------------------
    # LEIDEN
    # --------------------------------------------------

    print("\nRunning Leiden clustering...")

    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden",
        random_state=42,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    n_clusters = (
        adata.obs["leiden"]
        .nunique()
    )

    print(
        f"Leiden clusters: {n_clusters}"
    )

    # --------------------------------------------------
    # UMAP
    # --------------------------------------------------

    print("\nRunning UMAP...")

    sc.tl.umap(
        adata,
        random_state=42,
    )

    if "X_umap" not in adata.obsm:
        raise ValueError(
            "UMAP embedding was not generated"
        )

    # --------------------------------------------------
    # CLUSTER SUMMARY
    # --------------------------------------------------

    cluster_summary = (
        adata.obs["leiden"]
        .value_counts()
        .sort_index()
        .rename_axis("cluster")
        .reset_index(name="cells")
    )

    cluster_summary[
        "percentage"
    ] = (
        cluster_summary["cells"]
        / adata.n_obs
        * 100
    ).round(2)

    cluster_summary.to_csv(
        args.cluster_summary,
        index=False,
    )

    # --------------------------------------------------
    # SAMPLE × CLUSTER COMPOSITION
    # --------------------------------------------------

    sample_cluster = pd.crosstab(
        adata.obs["geo_accession"],
        adata.obs["leiden"],
    )

    sample_cluster.index.name = "sample"

    sample_cluster.to_csv(
        args.sample_cluster_summary
    )

    # --------------------------------------------------
    # UMAP BY SAMPLE
    # --------------------------------------------------

    sc.pl.umap(
        adata,
        color="geo_accession",
        show=False,
        title="CellFlowX — UMAP by Sample",
    )

    plt.savefig(
        args.umap_sample,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close()

    # --------------------------------------------------
    # UMAP BY CLUSTER
    # --------------------------------------------------

    sc.pl.umap(
        adata,
        color="leiden",
        show=False,
        legend_loc="on data",
        title="CellFlowX — Leiden Clusters",
    )

    plt.savefig(
        args.umap_cluster,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close()

    # --------------------------------------------------
    # METADATA
    # --------------------------------------------------

    adata.uns["cellflowx_embedding"] = {
        "pca_components": 50,
        "pca_hvg_only": True,
        "neighbors_n_neighbors":
            args.n_neighbors,
        "neighbors_n_pcs":
            args.n_pcs,
        "leiden_resolution":
            args.resolution,
        "random_state": 42,
        "integration_applied": False,
    }

    # --------------------------------------------------
    # SAVE
    # --------------------------------------------------

    output = Path(args.output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    adata.write_h5ad(
        output,
        compression="gzip",
    )

    # --------------------------------------------------
    # FINAL VALIDATION
    # --------------------------------------------------

    print(
        "\n=== PHASE 6 COMPLETE ==="
    )

    print(
        f"Cells: {adata.n_obs:,}"
    )

    print(
        f"Genes: {adata.n_vars:,}"
    )

    print(
        f"HVGs: {n_hvg:,}"
    )

    print(
        f"PCA dimensions: "
        f"{adata.obsm['X_pca'].shape}"
    )

    print(
        f"UMAP dimensions: "
        f"{adata.obsm['X_umap'].shape}"
    )

    print(
        f"Leiden clusters: "
        f"{n_clusters}"
    )

    print(
        f"Raw counts preserved: "
        f"{'counts' in adata.layers}"
    )

    print(
        f"Counts sparse: "
        f"{sparse.issparse(adata.layers['counts'])}"
    )

    print(
        f"\nSaved: {output.resolve()}"
    )


if __name__ == "__main__":
    main()
