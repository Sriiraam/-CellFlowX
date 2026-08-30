from pathlib import Path

import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data" / "processed" / "cellflowx_qc_metrics.h5ad"
OUTPUT = ROOT / "results" / "qc" / "qc_threshold_evaluation.csv"


THRESHOLDS = {
    "GSM8848584": {
        "min_genes": 500,
        "max_genes": 9500,
        "max_mt": 20,
    },
    "GSM8848585": {
        "min_genes": 500,
        "max_genes": 9000,
        "max_mt": 8,
    },
    "GSM8848586": {
        "min_genes": 450,
        "max_genes": 6500,
        "max_mt": 8,
    },
}


def main():

    adata = sc.read_h5ad(INPUT)

    records = []

    for sample, thresholds in THRESHOLDS.items():

        obs = adata.obs[
            adata.obs["geo_accession"] == sample
        ].copy()

        low_genes = (
            obs["n_genes_by_counts"] < thresholds["min_genes"]
        )

        high_genes = (
            obs["n_genes_by_counts"] > thresholds["max_genes"]
        )

        high_mt = (
            obs["pct_counts_mt"] > thresholds["max_mt"]
        )

        fail_any = low_genes | high_genes | high_mt

        total = len(obs)
        removed = int(fail_any.sum())
        retained = total - removed

        records.append({
            "sample": sample,
            "total_cells": total,

            "low_gene_cells": int(low_genes.sum()),
            "high_gene_cells": int(high_genes.sum()),
            "high_mt_cells": int(high_mt.sum()),

            "removed_any": removed,
            "retained_cells": retained,
            "retained_pct": round(retained / total * 100, 2),

            "min_genes": thresholds["min_genes"],
            "max_genes": thresholds["max_genes"],
            "max_mt": thresholds["max_mt"],
        })

    result = pd.DataFrame(records)

    result.to_csv(
        OUTPUT,
        index=False,
    )

    print("\n=== QC THRESHOLD EVALUATION ===\n")

    print(
        result[
            [
                "sample",
                "total_cells",
                "low_gene_cells",
                "high_gene_cells",
                "high_mt_cells",
                "removed_any",
                "retained_cells",
                "retained_pct",
            ]
        ].to_string(index=False)
    )

    print("\nTotal cells:", result["total_cells"].sum())
    print("Would remove:", result["removed_any"].sum())
    print("Would retain:", result["retained_cells"].sum())


if __name__ == "__main__":
    main()
