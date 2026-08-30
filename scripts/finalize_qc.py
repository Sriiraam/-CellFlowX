from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]

INPUT = (
    ROOT
    / "data"
    / "processed"
    / "cellflowx_qc_doublets.h5ad"
)

OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "cellflowx_qc.h5ad"
)

FLAGGED_OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "cellflowx_qc_flagged.h5ad"
)

SUMMARY_OUTPUT = (
    ROOT
    / "results"
    / "qc"
    / "final_qc_summary.csv"
)

REMOVAL_OUTPUT = (
    ROOT
    / "results"
    / "qc"
    / "qc_removal_reasons.csv"
)

DOUBLET_OUTPUT = (
    ROOT
    / "results"
    / "qc"
    / "doublet_candidates.csv"
)


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

    print(f"Loading: {INPUT}")

    adata = sc.read_h5ad(INPUT)

    # --------------------------------------------------
    # Basic validation
    # --------------------------------------------------

    assert adata.X is not None
    assert sparse.issparse(adata.X)
    assert adata.obs_names.is_unique
    assert adata.var_names.is_unique

    required = [
        "geo_accession",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "doublet_score",
        "predicted_doublet",
    ]

    missing = [
        col for col in required
        if col not in adata.obs.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required columns: {missing}"
        )

    # --------------------------------------------------
    # Initialize QC flags
    # --------------------------------------------------

    adata.obs["fail_low_genes"] = False
    adata.obs["fail_high_genes"] = False
    adata.obs["fail_high_mt"] = False

    # --------------------------------------------------
    # Apply sample-aware thresholds
    # --------------------------------------------------

    for sample, t in THRESHOLDS.items():

        mask = (
            adata.obs["geo_accession"] == sample
        )

        adata.obs.loc[
            mask,
            "fail_low_genes"
        ] = (
            adata.obs.loc[
                mask,
                "n_genes_by_counts"
            ] < t["min_genes"]
        )

        adata.obs.loc[
            mask,
            "fail_high_genes"
        ] = (
            adata.obs.loc[
                mask,
                "n_genes_by_counts"
            ] > t["max_genes"]
        )

        adata.obs.loc[
            mask,
            "fail_high_mt"
        ] = (
            adata.obs.loc[
                mask,
                "pct_counts_mt"
            ] > t["max_mt"]
        )

    # --------------------------------------------------
    # QC failure
    # --------------------------------------------------

    adata.obs["fail_basic_qc"] = (
        adata.obs["fail_low_genes"]
        | adata.obs["fail_high_genes"]
        | adata.obs["fail_high_mt"]
    )

    # --------------------------------------------------
    # Doublet handling
    #
    # Scrublet automatic global-rate estimates were
    # unstable, but only 7 cells were predicted.
    # We conservatively exclude those predicted cells.
    # --------------------------------------------------

    adata.obs["fail_doublet"] = (
        adata.obs["predicted_doublet"].astype(bool)
    )

    # --------------------------------------------------
    # Final QC decision
    # --------------------------------------------------

    adata.obs["qc_pass"] = ~(
        adata.obs["fail_basic_qc"]
        | adata.obs["fail_doublet"]
    )

    # --------------------------------------------------
    # Save doublet candidate table
    # --------------------------------------------------

    doublets = adata.obs[
        adata.obs["predicted_doublet"]
    ][
        [
            "geo_accession",
            "sample_id",
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "doublet_score",
            "fail_basic_qc",
        ]
    ].copy()

    doublets.to_csv(
        DOUBLET_OUTPUT
    )

    # --------------------------------------------------
    # Summary by sample
    # --------------------------------------------------

    records = []

    for sample in sorted(
        adata.obs["geo_accession"].unique()
    ):

        obs = adata.obs[
            adata.obs["geo_accession"] == sample
        ]

        total = len(obs)

        retained = int(
            obs["qc_pass"].sum()
        )

        removed = total - retained

        records.append({
            "sample": sample,
            "raw_cells": total,

            "low_gene_fail": int(
                obs["fail_low_genes"].sum()
            ),

            "high_gene_fail": int(
                obs["fail_high_genes"].sum()
            ),

            "high_mt_fail": int(
                obs["fail_high_mt"].sum()
            ),

            "predicted_doublets": int(
                obs["fail_doublet"].sum()
            ),

            "removed_total": removed,

            "retained_cells": retained,

            "retained_pct": round(
                retained / total * 100,
                2,
            ),
        })

    summary = pd.DataFrame(records)

    summary.to_csv(
        SUMMARY_OUTPUT,
        index=False,
    )

    # --------------------------------------------------
    # Removal reason table
    # --------------------------------------------------

    removed_obs = adata.obs[
        ~adata.obs["qc_pass"]
    ].copy()

    removal_cols = [
        "geo_accession",
        "sample_id",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "doublet_score",
        "predicted_doublet",
        "fail_low_genes",
        "fail_high_genes",
        "fail_high_mt",
        "fail_basic_qc",
        "fail_doublet",
    ]

    removed_obs[
        removal_cols
    ].to_csv(
        REMOVAL_OUTPUT
    )

    # --------------------------------------------------
    # Save flagged full object
    # --------------------------------------------------

    adata.write_h5ad(
        FLAGGED_OUTPUT,
        compression="gzip",
    )

    # --------------------------------------------------
    # Create final filtered object
    # --------------------------------------------------

    filtered = adata[
        adata.obs["qc_pass"]
    ].copy()

    # Remove unused QC-failure cells only.
    # Genes are deliberately NOT filtered here yet.
    # Gene filtering / preprocessing happens Phase 5.

    # --------------------------------------------------
    # Final validation
    # --------------------------------------------------

    assert filtered.X is not None
    assert sparse.issparse(filtered.X)
    assert filtered.obs_names.is_unique
    assert filtered.var_names.is_unique

    if filtered.n_obs == 0:
        raise ValueError(
            "QC removed all cells."
        )

    cell_counts = (
        filtered.obs["geo_accession"]
        .value_counts()
    )

    if (cell_counts == 0).any():
        raise ValueError(
            "One sample lost all cells."
        )

    X = filtered.X.tocsr()

    if X.nnz == 0:
        raise ValueError(
            "Filtered matrix contains no counts."
        )

    values = X.data

    integer_counts = np.allclose(
        values,
        np.round(values),
    )

    if not integer_counts:
        raise ValueError(
            "Count matrix is no longer integer-valued."
        )

    filtered.uns["qc_thresholds"] = THRESHOLDS

    filtered.uns["qc_notes"] = {
        "ribosomal_qc": (
            "Not used because canonical RPL/RPS "
            "genes are absent from supplied feature matrix."
        ),

        "doublet_detection": (
            "Scrublet applied independently per sample. "
            "Automatic overall-rate estimates were unstable; "
            "predicted doublet labels were conservatively "
            "excluded."
        ),

        "mitochondrial_genes_available": 11,

        "threshold_strategy": (
            "Sample-aware thresholds selected after "
            "inspection of sample-specific QC distributions "
            "and quantiles."
        ),
    }

    filtered.write_h5ad(
        OUTPUT,
        compression="gzip",
    )

    # --------------------------------------------------
    # Report
    # --------------------------------------------------

    print("\n=== FINAL QC SUMMARY ===\n")

    print(
        summary.to_string(
            index=False
        )
    )

    print("\n=== PROJECT TOTAL ===")

    raw = adata.n_obs
    final = filtered.n_obs
    removed = raw - final

    print(f"Raw cells:       {raw:,}")
    print(f"Removed cells:   {removed:,}")
    print(f"Retained cells:  {final:,}")
    print(
        f"Retention:       "
        f"{final / raw * 100:.2f}%"
    )

    print(
        f"Genes retained:  "
        f"{filtered.n_vars:,}"
    )

    print(
        f"Matrix sparse:   "
        f"{sparse.issparse(filtered.X)}"
    )

    print(
        f"Integer counts:  "
        f"{integer_counts}"
    )

    print("\nCells per sample:")

    print(
        filtered.obs[
            "geo_accession"
        ].value_counts().sort_index()
    )

    print("\nSaved:")
    print(OUTPUT)
    print(FLAGGED_OUTPUT)
    print(SUMMARY_OUTPUT)
    print(REMOVAL_OUTPUT)
    print(DOUBLET_OUTPUT)


if __name__ == "__main__":
    main()
