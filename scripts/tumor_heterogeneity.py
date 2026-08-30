import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--percentages", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--stacked-plot", required=True)
    p.add_argument("--heatmap", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)

    required = {"geo_accession", "cell_type"}

    missing = required - set(adata.obs.columns)

    if missing:
        raise ValueError(
            f"Missing metadata: {sorted(missing)}"
        )

    print(
        f"Input: {adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    # --------------------------------------------
    # CELL COUNTS
    # --------------------------------------------

    counts = pd.crosstab(
        adata.obs["geo_accession"],
        adata.obs["cell_type"]
    )

    counts.index.name = "sample"

    counts.to_csv(args.counts)

    # --------------------------------------------
    # WITHIN-SAMPLE PERCENTAGES
    # --------------------------------------------

    percentages = (
        counts.div(
            counts.sum(axis=1),
            axis=0
        ) * 100
    )

    percentages = percentages.round(2)

    percentages.to_csv(args.percentages)

    # --------------------------------------------
    # SUMMARY
    # --------------------------------------------

    records = []

    for sample in percentages.index:

        row = percentages.loc[sample]

        dominant_state = row.idxmax()
        dominant_pct = row.max()

        detected_states = int(
            (counts.loc[sample] > 0).sum()
        )

        records.append({
            "sample": sample,
            "total_cells":
                int(counts.loc[sample].sum()),
            "detected_states":
                detected_states,
            "dominant_state":
                dominant_state,
            "dominant_state_pct":
                round(float(dominant_pct), 2)
        })

    summary = pd.DataFrame(records)

    summary.to_csv(
        args.summary,
        index=False
    )

    # --------------------------------------------
    # STACKED COMPOSITION PLOT
    # --------------------------------------------

    ax = percentages.plot(
        kind="bar",
        stacked=True,
        figsize=(11, 6)
    )

    ax.set_ylabel("Cells (%)")
    ax.set_xlabel("Metastatic sample")
    ax.set_title(
        "CellFlowX — Inter-Tumor Cellular Composition"
    )

    ax.legend(
        title="Cell type / state",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    plt.xticks(
        rotation=0
    )

    plt.tight_layout()

    plt.savefig(
        args.stacked_plot,
        dpi=220,
        bbox_inches="tight"
    )

    plt.close()

    # --------------------------------------------
    # HEATMAP
    # --------------------------------------------

    fig, ax = plt.subplots(
        figsize=(11, 4)
    )

    matrix = percentages.to_numpy()

    image = ax.imshow(
        matrix,
        aspect="auto"
    )

    ax.set_xticks(
        range(len(percentages.columns))
    )

    ax.set_xticklabels(
        percentages.columns,
        rotation=60,
        ha="right"
    )

    ax.set_yticks(
        range(len(percentages.index))
    )

    ax.set_yticklabels(
        percentages.index
    )

    ax.set_title(
        "CellFlowX — Cell-State Composition Heatmap"
    )

    ax.set_xlabel(
        "Cell type / transcriptional state"
    )

    ax.set_ylabel(
        "Metastatic sample"
    )

    fig.colorbar(
        image,
        ax=ax,
        label="Cells (%)"
    )

    fig.tight_layout()

    fig.savefig(
        args.heatmap,
        dpi=220,
        bbox_inches="tight"
    )

    plt.close(fig)

    # --------------------------------------------
    # PRINT RESULTS
    # --------------------------------------------

    print("\n=== CELL-TYPE PERCENTAGES ===")
    print(percentages.to_string())

    print("\n=== HETEROGENEITY SUMMARY ===")
    print(summary.to_string(index=False))

    print("\n=== PHASE 8 COMPLETE ===")

    for path in [
        args.counts,
        args.percentages,
        args.summary,
        args.stacked_plot,
        args.heatmap
    ]:
        print(Path(path).resolve())


if __name__ == "__main__":
    main()
