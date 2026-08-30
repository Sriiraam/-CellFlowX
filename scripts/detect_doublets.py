from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr


ROOT = Path(__file__).resolve().parents[1]

INPUT = ROOT / "data" / "processed" / "cellflowx_qc_metrics.h5ad"

OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "cellflowx_qc_doublets.h5ad"
)

SUMMARY = (
    ROOT
    / "results"
    / "qc"
    / "doublet_summary.csv"
)


def main():

    adata = sc.read_h5ad(INPUT)

    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False

    records = []

    for sample in sorted(
        adata.obs["geo_accession"].unique()
    ):

        print(f"\n===== {sample} =====")

        mask = (
            adata.obs["geo_accession"] == sample
        )

        subset = adata[mask].copy()

        scrub = scr.Scrublet(
            subset.X,
            expected_doublet_rate=0.06,
            random_state=42,
        )

        scores, predicted = scrub.scrub_doublets()

        adata.obs.loc[
            mask,
            "doublet_score"
        ] = scores

        adata.obs.loc[
            mask,
            "predicted_doublet"
        ] = predicted

        n_cells = subset.n_obs
        n_doublets = int(predicted.sum())

        print("Cells:", n_cells)
        print("Predicted doublets:", n_doublets)
        print(
            "Predicted doublet %:",
            round(n_doublets / n_cells * 100, 2),
        )
        print("Threshold:", scrub.threshold_)

        records.append({
            "sample": sample,
            "cells": n_cells,
            "predicted_doublets": n_doublets,
            "predicted_doublet_pct": round(
                n_doublets / n_cells * 100,
                2,
            ),
            "scrublet_threshold": scrub.threshold_,
        })

    summary = pd.DataFrame(records)

    print("\n=== DOUBLET SUMMARY ===")
    print(summary.to_string(index=False))

    SUMMARY.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary.to_csv(
        SUMMARY,
        index=False,
    )

    adata.write_h5ad(
        OUTPUT,
        compression="gzip",
    )

    print(f"\nSaved: {OUTPUT}")
    print(f"Saved: {SUMMARY}")


if __name__ == "__main__":
    main()
