import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)

    with open(args.thresholds) as f:
        thresholds_all = json.load(f)

    records = []

    for sample, thresholds in thresholds_all.items():

        obs = adata.obs[
            adata.obs["geo_accession"] == sample
        ]

        low_genes = (
            obs["n_genes_by_counts"]
            < thresholds["min_genes"]
        )

        high_genes = (
            obs["n_genes_by_counts"]
            > thresholds["max_genes"]
        )

        high_mt = (
            obs["pct_counts_mt"]
            > thresholds["max_mt"]
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
            "retained_pct": round(
                retained / total * 100,
                2
            ),
        })

    result = pd.DataFrame(records)

    output = Path(args.output)

    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result.to_csv(
        output,
        index=False
    )

    print("\n=== QC THRESHOLD EVALUATION ===\n")
    print(result.to_string(index=False))

    print("\nTotal cells:", result["total_cells"].sum())
    print("Would remove:", result["removed_any"].sum())
    print("Would retain:", result["retained_cells"].sum())

    print(f"\nSaved: {output.resolve()}")


if __name__ == "__main__":
    main()
