import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import gseapy as gp


STATES = [
    "Prostate epithelial - AR high",
    "Prostate epithelial - luminal",
    "Neuroendocrine-like",
    "Steroidogenic-like",
]


COMPARISONS = [
    (
        "Prostate epithelial - AR high",
        "Neuroendocrine-like",
        "AR_high_vs_Neuroendocrine"
    ),
    (
        "Prostate epithelial - luminal",
        "Neuroendocrine-like",
        "Luminal_vs_Neuroendocrine"
    ),
    (
        "Prostate epithelial - AR high",
        "Steroidogenic-like",
        "AR_high_vs_Steroidogenic"
    ),
    (
        "Neuroendocrine-like",
        "Steroidogenic-like",
        "Neuroendocrine_vs_Steroidogenic"
    ),
]


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def safe_symbol_map(adata):
    if "gene_symbols" not in adata.var.columns:
        raise ValueError("adata.var['gene_symbols'] missing")

    return pd.Series(
        adata.var["gene_symbols"].astype(str).values,
        index=adata.var_names.astype(str)
    ).to_dict()


def volcano(df, title, output):
    x = df["logfoldchanges"].replace(
        [np.inf, -np.inf], np.nan
    )

    p = df["pvals_adj"].clip(lower=1e-300)
    y = -np.log10(p)

    significant = (
        (df["pvals_adj"] < 0.05) &
        (df["logfoldchanges"].abs() >= 1)
    )

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        x[~significant],
        y[~significant],
        s=7,
        alpha=0.35
    )

    ax.scatter(
        x[significant],
        y[significant],
        s=9,
        alpha=0.7
    )

    ax.axvline(1, linestyle="--", linewidth=1)
    ax.axvline(-1, linestyle="--", linewidth=1)
    ax.axhline(
        -np.log10(0.05),
        linestyle="--",
        linewidth=1
    )

    ax.set_xlabel("log fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(title)

    # Label strongest genes
    label_df = (
        df.loc[significant]
        .assign(abs_lfc=lambda z: z["logfoldchanges"].abs())
        .sort_values(
            ["pvals_adj", "abs_lfc"],
            ascending=[True, False]
        )
        .head(10)
    )

    for _, row in label_df.iterrows():
        if pd.notna(row["gene_symbol"]):
            ax.annotate(
                row["gene_symbol"],
                (
                    row["logfoldchanges"],
                    -np.log10(max(row["pvals_adj"], 1e-300))
                ),
                fontsize=7
            )

    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def enrichment(de, label, output_dir):
    genes = (
        de.loc[
            (de["pvals_adj"] < 0.05) &
            (de["logfoldchanges"] >= 1),
            "gene_symbol"
        ]
        .dropna()
        .astype(str)
    )

    genes = genes[
        ~genes.isin(["nan", "None", ""])
    ].drop_duplicates()

    if len(genes) < 10:
        print(
            f"{label}: only {len(genes)} upregulated genes; "
            "skipping enrichment"
        )
        return

    print(
        f"{label}: enrichment using "
        f"{len(genes)} upregulated genes"
    )

    try:
        enr = gp.enrichr(
            gene_list=genes.tolist(),
            gene_sets="GO_Biological_Process_2023",
            organism="Human",
            outdir=None,
            cutoff=0.05
        )

        result = enr.results.copy()

        if result.empty:
            print(f"{label}: no significant pathways")
            return

        result.to_csv(
            output_dir / f"{label}_GO_BP.csv",
            index=False
        )

        top = (
            result
            .sort_values("Adjusted P-value")
            .head(15)
            .copy()
        )

        top["minus_log10_fdr"] = (
            -np.log10(
                top["Adjusted P-value"]
                .clip(lower=1e-300)
            )
        )

        fig, ax = plt.subplots(figsize=(9, 6))

        ax.barh(
            top["Term"][::-1],
            top["minus_log10_fdr"][::-1]
        )

        ax.set_xlabel("-log10 adjusted p-value")
        ax.set_title(
            f"{label} — GO Biological Process"
        )

        fig.tight_layout()

        fig.savefig(
            output_dir / f"{label}_GO_BP.png",
            dpi=220,
            bbox_inches="tight"
        )

        plt.close(fig)

    except Exception as e:
        # DE remains valid even if internet-based
        # enrichment service is unavailable.
        print(
            f"WARNING: enrichment failed for {label}: {e}"
        )


def main():
    a = args()

    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(a.input)

    print(
        f"Input: {adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    if "cell_type" not in adata.obs.columns:
        raise ValueError("cell_type missing")

    available = set(
        adata.obs["cell_type"].astype(str).unique()
    )

    missing = set(STATES) - available

    if missing:
        raise ValueError(
            f"Missing tumor states: {sorted(missing)}"
        )

    symbol_map = safe_symbol_map(adata)

    # -------------------------------------------------
    # Sample/state representation
    # -------------------------------------------------

    representation = pd.crosstab(
        adata.obs["geo_accession"],
        adata.obs["cell_type"]
    )

    representation = representation.reindex(
        columns=STATES,
        fill_value=0
    )

    representation.to_csv(
        out / "tumor_state_sample_representation.csv"
    )

    print("\n=== SAMPLE × TUMOR STATE ===")
    print(representation.to_string())

    # -------------------------------------------------
    # Restrict to tumor-associated states
    # -------------------------------------------------

    mask = (
        adata.obs["cell_type"]
        .astype(str)
        .isin(STATES)
    )

    tumor = adata[mask].copy()

    print(
        "\nTumor-associated/state cells:",
        f"{tumor.n_obs:,}"
    )

    # -------------------------------------------------
    # Pairwise exploratory DE
    # -------------------------------------------------

    summary_records = []

    for state_a, state_b, label in COMPARISONS:

        print(f"\n=== {label} ===")

        pair_mask = (
            tumor.obs["cell_type"]
            .astype(str)
            .isin([state_a, state_b])
        )

        pair = tumor[pair_mask].copy()

        n_a = int(
            (
                pair.obs["cell_type"].astype(str)
                == state_a
            ).sum()
        )

        n_b = int(
            (
                pair.obs["cell_type"].astype(str)
                == state_b
            ).sum()
        )

        print(state_a, n_a)
        print(state_b, n_b)

        if n_a < 20 or n_b < 20:
            print("Insufficient cells — skipping")
            continue

        sc.tl.rank_genes_groups(
            pair,
            groupby="cell_type",
            groups=[state_a],
            reference=state_b,
            method="wilcoxon",
            use_raw=False,
            pts=True
        )

        de = sc.get.rank_genes_groups_df(
            pair,
            group=state_a
        )

        de["gene_symbol"] = (
            de["names"]
            .astype(str)
            .map(symbol_map)
        )

        de["comparison"] = label
        de["state_A"] = state_a
        de["state_B"] = state_b

        de.to_csv(
            out / f"{label}_DE.csv",
            index=False
        )

        significant = de[
            (de["pvals_adj"] < 0.05) &
            (de["logfoldchanges"].abs() >= 1)
        ]

        up = int(
            (
                significant["logfoldchanges"] >= 1
            ).sum()
        )

        down = int(
            (
                significant["logfoldchanges"] <= -1
            ).sum()
        )

        summary_records.append({
            "comparison": label,
            "state_A": state_a,
            "state_B": state_b,
            "cells_A": n_a,
            "cells_B": n_b,
            "significant_genes": len(significant),
            "up_in_A": up,
            "up_in_B": down
        })

        volcano(
            de,
            label.replace("_", " "),
            out / f"{label}_volcano.png"
        )

        # Pathways enriched among genes higher in A
        enrichment(
            de,
            f"{label}_up_A",
            out
        )

        # Reverse sign to examine genes higher in B
        reverse = de.copy()
        reverse["logfoldchanges"] *= -1

        enrichment(
            reverse,
            f"{label}_up_B",
            out
        )

    summary = pd.DataFrame(summary_records)

    summary.to_csv(
        out / "tumor_state_DE_summary.csv",
        index=False
    )

    print("\n=== DE SUMMARY ===")
    print(summary.to_string(index=False))

    print("\n=== PHASE 10 COMPLETE ===")


if __name__ == "__main__":
    main()
