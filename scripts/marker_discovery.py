import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


MARKER_PANEL = {
    "Epithelial": [
        "EPCAM", "KRT8", "KRT18", "KRT19"
    ],
    "Prostate_epithelial": [
        "AR", "KLK3", "KLK2", "ACPP", "NKX3-1"
    ],
    "T_cell": [
        "CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2"
    ],
    "CD8_T": [
        "CD8A", "CD8B"
    ],
    "CD4_T": [
        "CD4", "IL7R", "LTB"
    ],
    "NK": [
        "NKG7", "GNLY", "KLRD1", "PRF1"
    ],
    "B_cell": [
        "CD79A", "MS4A1", "CD37", "CD74"
    ],
    "Plasma": [
        "JCHAIN", "MZB1", "SDC1"
    ],
    "Myeloid": [
        "LYZ", "CTSS", "FCER1G", "TYROBP"
    ],
    "Macrophage": [
        "C1QA", "C1QB", "C1QC", "APOE"
    ],
    "Endothelial": [
        "PECAM1", "VWF", "EMCN", "KDR"
    ],
    "Fibroblast": [
        "COL1A1", "COL1A2", "DCN", "COL3A1"
    ]
}


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--input", required=True)
    p.add_argument("--markers", required=True)
    p.add_argument("--top-markers", required=True)
    p.add_argument("--dotplot", required=True)

    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)

    if "leiden" not in adata.obs:
        raise ValueError("Leiden clusters missing")

    print(
        f"Input: {adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    # --------------------------------------------------
    # Differential marker discovery
    # --------------------------------------------------

    print("\nFinding cluster markers...")

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        use_raw=False,
        pts=True,
    )

    markers = sc.get.rank_genes_groups_df(
        adata,
        group=None,
    )

    markers.to_csv(
        args.markers,
        index=False,
    )

    # Stronger shortlist
    top = markers[
        (markers["logfoldchanges"] > 0.5) &
        (markers["pvals_adj"] < 0.05)
    ].copy()

    top = (
        top.sort_values(
            ["group", "scores"],
            ascending=[True, False]
        )
        .groupby("group", observed=True)
        .head(20)
    )

    top.to_csv(
        args.top_markers,
        index=False,
    )

    # --------------------------------------------------
    # Canonical marker panel
    # --------------------------------------------------

    symbols = set(
        adata.var["gene_symbols"]
        .astype(str)
    )

    available = {}

    print("\n=== CANONICAL MARKERS AVAILABLE ===")

    for cell_type, genes in MARKER_PANEL.items():

        present = [
            gene for gene in genes
            if gene in symbols
        ]

        available[cell_type] = present

        print(
            f"{cell_type}: "
            f"{', '.join(present) if present else 'NONE'}"
        )

    # Scanpy plotting uses var_names.
    # Current var_names are Ensembl IDs, so map symbols → IDs.

    symbol_to_id = {}

    for gene_id, symbol in zip(
        adata.var_names,
        adata.var["gene_symbols"].astype(str)
    ):
        if symbol not in symbol_to_id:
            symbol_to_id[symbol] = gene_id

    plot_markers = {}

    for cell_type, genes in available.items():

        ids = [
            symbol_to_id[g]
            for g in genes
            if g in symbol_to_id
        ]

        if ids:
            plot_markers[cell_type] = ids

    # --------------------------------------------------
    # Dotplot
    # --------------------------------------------------

    sc.pl.dotplot(
        adata,
        var_names=plot_markers,
        groupby="leiden",
        standard_scale="var",
        show=False,
        save=None,
    )

    import matplotlib.pyplot as plt

    plt.savefig(
        args.dotplot,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close("all")

    print("\n=== MARKER DISCOVERY COMPLETE ===")
    print(
        "Clusters:",
        adata.obs["leiden"].nunique()
    )

    print(
        f"Marker rows: {len(markers):,}"
    )

    print(f"Saved: {Path(args.markers).resolve()}")
    print(f"Saved: {Path(args.top_markers).resolve()}")
    print(f"Saved: {Path(args.dotplot).resolve()}")


if __name__ == "__main__":
    main()
