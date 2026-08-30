import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)

    return parser.parse_args()


def main():
    args = parse_args()

    print(f"Loading: {args.input}")

    adata = sc.read_h5ad(args.input)

    if adata.X is None:
        raise ValueError("adata.X is missing")

    if not sparse.issparse(adata.X):
        raise ValueError("Expected sparse matrix")

    X = adata.X.tocsr()

    symbols = (
        adata.var["gene_symbols"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    mt_mask = symbols.str.startswith("MT-").to_numpy()

    print("\nQC gene sets:")
    print(f"Mitochondrial genes: {mt_mask.sum()}")
    print(f"RPL genes: {symbols.str.startswith('RPL').sum()}")
    print(f"RPS genes: {symbols.str.startswith('RPS').sum()}")

    total_counts = np.asarray(
        X.sum(axis=1)
    ).ravel()

    n_genes = np.diff(X.indptr)

    mt_counts = np.asarray(
        X[:, mt_mask].sum(axis=1)
    ).ravel()

    pct_mt = np.divide(
        mt_counts,
        total_counts,
        out=np.zeros_like(mt_counts, dtype=float),
        where=total_counts > 0,
    ) * 100

    adata.obs["total_counts"] = total_counts
    adata.obs["n_genes_by_counts"] = n_genes
    adata.obs["total_counts_mt"] = mt_counts
    adata.obs["pct_counts_mt"] = pct_mt

    metrics = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
    ]

    print("\n=== GLOBAL QC SUMMARY ===")
    print(
        adata.obs[metrics]
        .describe()
        .round(2)
    )

    summary = (
        adata.obs
        .groupby(
            "geo_accession",
            observed=True,
        )
        .agg(
            cells=("geo_accession", "size"),
            median_counts=("total_counts", "median"),
            median_genes=("n_genes_by_counts", "median"),
            median_pct_mt=("pct_counts_mt", "median"),
        )
        .round(2)
    )

    print("\n=== QC BY SAMPLE ===")
    print(summary)

    output = Path(args.output)
    summary_path = Path(args.summary)

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary.to_csv(summary_path)

    adata.uns["qc_notes"] = {
        "mitochondrial_genes_detected":
            int(mt_mask.sum()),

        "ribosomal_qc_available":
            False,

        "ribosomal_qc_reason":
            "Canonical RPL/RPS genes absent from supplied matrix",
    }

    adata.write_h5ad(
        output,
        compression="gzip",
    )

    print(f"\nSaved AnnData: {output.resolve()}")
    print(f"Saved summary: {summary_path.resolve()}")


if __name__ == "__main__":
    main()
