import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--umap", required=True)
    p.add_argument("--cluster-table", required=True)
    p.add_argument("--sample-table", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)

    if "leiden" not in adata.obs:
        raise ValueError("Leiden clusters missing")

    with open(args.annotations) as f:
        mapping = json.load(f)

    clusters = set(
        adata.obs["leiden"].astype(str).unique()
    )

    missing = clusters - set(mapping)

    if missing:
        raise ValueError(
            f"Missing annotations for clusters: {sorted(missing)}"
        )

    adata.obs["cell_type"] = (
        adata.obs["leiden"]
        .astype(str)
        .map(mapping)
        .astype("category")
    )

    # Cluster annotation table
    cluster_table = (
        adata.obs[
            ["leiden", "cell_type"]
        ]
        .drop_duplicates()
        .sort_values(
            "leiden",
            key=lambda x: x.astype(int)
        )
    )

    cluster_table.to_csv(
        args.cluster_table,
        index=False
    )

    # Sample × cell-type composition
    sample_table = pd.crosstab(
        adata.obs["geo_accession"],
        adata.obs["cell_type"]
    )

    sample_table.index.name = "sample"

    sample_table.to_csv(
        args.sample_table
    )

    # Annotated UMAP
    sc.pl.umap(
        adata,
        color="cell_type",
        show=False,
        legend_loc="right margin",
        title="CellFlowX — Cell-Type Annotation"
    )

    plt.savefig(
        args.umap,
        dpi=220,
        bbox_inches="tight"
    )

    plt.close("all")

    adata.uns["annotation"] = {
        "method":
            "manual evidence-based cluster annotation",
        "basis":
            "cluster differential markers and canonical marker panel",
        "malignancy_status":
            "not definitively assigned; CNV evidence not yet evaluated"
    }

    output = Path(args.output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    adata.write_h5ad(
        output,
        compression="gzip"
    )

    print("\n=== PHASE 7 ANNOTATION COMPLETE ===")
    print("Cells:", f"{adata.n_obs:,}")
    print("Clusters:", adata.obs["leiden"].nunique())
    print("Annotations:", adata.obs["cell_type"].nunique())

    print("\nCell types:")
    print(
        adata.obs["cell_type"]
        .value_counts()
        .to_string()
    )

    print(f"\nSaved: {output.resolve()}")


if __name__ == "__main__":
    main()
