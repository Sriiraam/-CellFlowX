import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]






def save_histogram(adata, metric, xlabel, filename):
    fig, ax = plt.subplots(figsize=(8, 5))

    for sample in sorted(adata.obs["geo_accession"].unique()):
        values = adata.obs.loc[
            adata.obs["geo_accession"] == sample,
            metric,
        ]

        ax.hist(
            values,
            bins=60,
            alpha=0.45,
            label=sample,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cells")
    ax.legend()

    fig.tight_layout()
    fig.savefig(
        OUTDIR / filename,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


def main():
    global OUTDIR

    parser = argparse.ArgumentParser(
        description="Generate CellFlowX QC plots."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    INPUT = Path(args.input)
    OUTDIR = Path(args.outdir)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(INPUT)

    metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]

    # -------------------------------
    # Quantiles by sample
    # -------------------------------

    records = []

    quantiles = [
        0.01,
        0.05,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.95,
        0.99,
    ]

    for sample in sorted(adata.obs["geo_accession"].unique()):

        subset = adata.obs[
            adata.obs["geo_accession"] == sample
        ]

        for metric in metrics:

            q = subset[metric].quantile(quantiles)

            for quantile, value in q.items():

                records.append(
                    {
                        "sample": sample,
                        "metric": metric,
                        "quantile": quantile,
                        "value": value,
                    }
                )

    quantile_df = pd.DataFrame(records)

    quantile_df.to_csv(
        OUTDIR / "qc_quantiles_by_sample.csv",
        index=False,
    )

    print("\n=== QC QUANTILES ===")

    for metric in metrics:

        print(f"\n--- {metric} ---")

        table = (
            quantile_df[
                quantile_df["metric"] == metric
            ]
            .pivot(
                index="sample",
                columns="quantile",
                values="value",
            )
            .round(2)
        )

        print(table)

    # -------------------------------
    # Histograms
    # -------------------------------

    save_histogram(
        adata,
        "total_counts",
        "Total counts per cell",
        "total_counts_distribution.png",
    )

    save_histogram(
        adata,
        "n_genes_by_counts",
        "Detected genes per cell",
        "genes_per_cell_distribution.png",
    )

    save_histogram(
        adata,
        "pct_counts_mt",
        "Mitochondrial counts (%)",
        "mitochondrial_distribution.png",
    )

    # -------------------------------
    # Counts vs genes scatter
    # -------------------------------

    fig, ax = plt.subplots(figsize=(8, 6))

    scatter = ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        s=8,
        alpha=0.5,
    )

    ax.set_xlabel("Total counts")
    ax.set_ylabel("Detected genes")

    fig.colorbar(
        scatter,
        ax=ax,
        label="Mitochondrial %",
    )

    fig.tight_layout()

    fig.savefig(
        OUTDIR / "counts_vs_genes.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    print("\nSaved QC plots to:", OUTDIR)


if __name__ == "__main__":
    main()
