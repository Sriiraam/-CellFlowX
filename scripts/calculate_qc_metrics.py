from pathlib import Path

import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "cellflowx_raw_merged.h5ad"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "cellflowx_qc_metrics.h5ad"
)

SUMMARY_FILE = (
    PROJECT_ROOT
    / "results"
    / "qc"
    / "qc_summary_by_sample.csv"
)


def main():

    print(f"Loading: {INPUT_FILE}")

    adata = sc.read_h5ad(INPUT_FILE)

    # --------------------------------
    # Gene annotations for QC
    # --------------------------------

    symbols = adata.var["gene_symbols"].astype(str)

    adata.var["mt"] = symbols.str.upper().str.startswith("MT-")
    adata.var["ribo"] = symbols.str.upper().str.startswith(
        ("RPS", "RPL")
    )

    print("\nQC gene sets:")
    print(f"Mitochondrial genes: {adata.var['mt'].sum()}")
    print(f"Ribosomal genes: {adata.var['ribo'].sum()}")

    # --------------------------------
    # Calculate QC metrics
    # --------------------------------

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # --------------------------------
    # Sanity checks
    # --------------------------------

    required_metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
    ]

    for metric in required_metrics:
        if metric not in adata.obs.columns:
            raise ValueError(f"Missing QC metric: {metric}")

    print("\n=== GLOBAL QC SUMMARY ===")

    print(
        adata.obs[
            required_metrics
        ].describe().round(2)
    )

    # --------------------------------
    # Per-sample summary
    # --------------------------------

    summary = (
        adata.obs
        .groupby("geo_accession", observed=True)
        .agg(
            cells=("geo_accession", "size"),

            median_counts=("total_counts", "median"),
            median_genes=("n_genes_by_counts", "median"),

            median_pct_mt=("pct_counts_mt", "median"),
            median_pct_ribo=("pct_counts_ribo", "median"),
        )
        .round(2)
    )

    print("\n=== QC BY SAMPLE ===")
    print(summary)

    # --------------------------------
    # Save
    # --------------------------------

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_FILE.parent.mkdir(parents=True, exist_ok=True)

    summary.to_csv(SUMMARY_FILE)

    adata.write_h5ad(
        OUTPUT_FILE,
        compression="gzip",
    )

    print(f"\nSaved AnnData: {OUTPUT_FILE}")
    print(f"Saved summary: {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
