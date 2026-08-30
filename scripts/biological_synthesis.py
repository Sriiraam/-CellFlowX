import argparse
from pathlib import Path
import shutil
import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--heterogeneity-dir", required=True)
    p.add_argument("--cnv-dir", required=True)
    p.add_argument("--enrichment-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    het_dir = Path(args.heterogeneity_dir)
    cnv_dir = Path(args.cnv_dir)
    enr_dir = Path(args.enrichment_dir)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # 1. SAMPLE-LEVEL BIOLOGICAL SUMMARY
    # --------------------------------------------------
    heterogeneity = pd.read_csv(
        het_dir / "heterogeneity_summary.csv"
    )

    cnv = pd.read_csv(
        cnv_dir / "cnv_summary.csv"
    )

    cnv_candidates = cnv[
        cnv["cell_type"].isin([
            "Prostate epithelial - AR high",
            "Prostate epithelial - luminal",
            "Neuroendocrine-like",
            "Steroidogenic-like",
            "Epithelial-like - uncertain",
        ])
    ].copy()

    idx = (
        cnv_candidates
        .groupby("geo_accession")["cnv_high_pct"]
        .idxmax()
    )

    strongest_cnv = cnv_candidates.loc[
        idx,
        [
            "geo_accession",
            "cell_type",
            "median_cnv_score",
            "cnv_high_pct"
        ]
    ].copy()

    strongest_cnv.columns = [
        "sample",
        "strongest_cnv_candidate_state",
        "median_cnv_score",
        "cnv_high_pct"
    ]

    biological_summary = heterogeneity.rename(
        columns={"sample": "sample"}
    ).merge(
        strongest_cnv,
        on="sample",
        how="left"
    )

    biological_summary.to_csv(
        out / "biological_summary.csv",
        index=False
    )

    # --------------------------------------------------
    # 2. PATHWAY SUMMARY
    # --------------------------------------------------
    enrichment = pd.read_csv(
        enr_dir / "functional_enrichment_summary.csv"
    )

    pathway_summary = enrichment[
        [
            "gene_set",
            "input_genes",
            "significant_pathways",
            "top_pathway",
            "top_adjusted_p"
        ]
    ].copy()

    pathway_summary.to_csv(
        out / "state_pathway_summary.csv",
        index=False
    )

    # --------------------------------------------------
    # 3. FINAL FLAGSHIP FIGURE
    # --------------------------------------------------
    composition = pd.read_csv(
        het_dir / "celltype_percentages_by_sample.csv"
    )

    # Phase 8 composition table uses "sample" as the sample identifier.
    # Accept geo_accession as well for compatibility.
    if "geo_accession" in composition.columns:
        composition = composition.set_index("geo_accession")
    elif "sample" in composition.columns:
        composition = composition.set_index("sample")
    else:
        raise KeyError(
            f"Could not identify sample column. Columns: {list(composition.columns)}"
        )

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.2, 1],
        height_ratios=[1, 1],
        hspace=0.35,
        wspace=0.30
    )

    # Panel A — composition
    ax1 = fig.add_subplot(gs[0, 0])

    composition.plot(
        kind="bar",
        stacked=True,
        ax=ax1,
        width=0.75
    )

    ax1.set_title(
        "A. Inter-tumor cellular composition",
        loc="left",
        fontweight="bold"
    )
    ax1.set_xlabel("")
    ax1.set_ylabel("Cell composition (%)")
    ax1.tick_params(axis="x", rotation=0)

    ax1.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=7,
        frameon=False
    )

    # Panel B — CNV-like evidence
    ax2 = fig.add_subplot(gs[0, 1])

    cnv_plot = cnv_candidates.pivot_table(
        index="cell_type",
        columns="geo_accession",
        values="cnv_high_pct",
        fill_value=0
    )

    im = ax2.imshow(
        cnv_plot.values,
        aspect="auto"
    )

    ax2.set_xticks(range(len(cnv_plot.columns)))
    ax2.set_xticklabels(
        cnv_plot.columns,
        rotation=30,
        ha="right"
    )

    ax2.set_yticks(range(len(cnv_plot.index)))
    ax2.set_yticklabels(
        cnv_plot.index,
        fontsize=8
    )

    ax2.set_title(
        "B. CNV-like expression evidence",
        loc="left",
        fontweight="bold"
    )

    cbar = fig.colorbar(im, ax=ax2)
    cbar.set_label("CNV-high cells (%)")

    # Panel C — dominant states
    ax3 = fig.add_subplot(gs[1, 0])

    dominant = biological_summary[
        ["sample", "dominant_state", "dominant_state_pct"]
    ].copy()

    bars = ax3.bar(
        dominant["sample"],
        dominant["dominant_state_pct"]
    )

    ax3.set_ylim(
        0,
        max(65, dominant["dominant_state_pct"].max() + 8)
    )

    ax3.set_ylabel("Dominant state (%)")
    ax3.set_title(
        "C. Dominant cellular state by tumor",
        loc="left",
        fontweight="bold"
    )

    for bar, state, pct in zip(
        bars,
        dominant["dominant_state"],
        dominant["dominant_state_pct"]
    ):
        ax3.text(
            bar.get_x() + bar.get_width()/2,
            bar.get_height() + 1,
            f"{state}\n{pct:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8
        )

    # Panel D — top pathway evidence
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")

    pathway_lines = []

    for _, row in pathway_summary.iterrows():
        comparison = row["gene_set"].replace("_", " ")
        pathway = row["top_pathway"]

        pathway_lines.append(
            f"• {comparison}\n  {pathway}"
        )

    pathway_text = "\n\n".join(pathway_lines)

    ax4.text(
        0,
        1,
        "D. Leading functional programs",
        fontsize=12,
        fontweight="bold",
        va="top"
    )

    ax4.text(
        0,
        0.91,
        pathway_text,
        fontsize=8,
        va="top",
        wrap=True
    )

    fig.suptitle(
        "CellFlowX — Biological Synthesis of Metastatic Prostate Cancer Heterogeneity",
        fontsize=16,
        fontweight="bold"
    )

    plt.savefig(
        out / "cellflowx_biological_synthesis.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # --------------------------------------------------
    # 4. KEY FINDINGS MARKDOWN
    # --------------------------------------------------
    findings = """# CellFlowX — Key Biological Findings

## Dataset

CellFlowX profiles three metastatic prostate cancer tumor samples using
single-cell RNA sequencing.

A total of 8,233 high-quality cells were retained after quality control.

## Inter-tumor heterogeneity

Strong differences in cellular composition were observed across the three
tumor samples.

- GSM8848584 was dominated by steroidogenic-like and AR-high prostate
  epithelial populations.
- GSM8848585 was dominated by neuroendocrine-like cells, with substantial
  luminal and cycling populations.
- GSM8848586 was dominated by activated fibroblasts and macrophages,
  indicating a strongly stromal/myeloid-rich microenvironment.

## CNV-like transcriptional evidence

Chromosome-level expression-deviation analysis identified stronger CNV-like
signals in selected prostate epithelial populations.

The strongest evidence included luminal prostate epithelial cells in
GSM8848585 and AR-high prostate epithelial cells in GSM8848584.

These measurements are expression-based CNV proxies and must not be
interpreted as DNA-level copy-number calls.

## Tumor-associated transcriptional states

Exploratory differential-expression analysis identified substantial
transcriptional separation among AR-high, luminal, neuroendocrine-like and
steroidogenic-like states.

Neuroendocrine-like populations showed enrichment for neuronal and synaptic
programs, including nervous-system development and chemical synaptic
transmission.

Steroidogenic-like populations showed strong enrichment for cholesterol
metabolic processes.

## Interpretation

Together, the results demonstrate pronounced inter-tumor heterogeneity in
cellular composition and transcriptional state across the metastatic prostate
cancer samples analyzed by CellFlowX.

## Limitations

Only three metastatic tumor samples were analyzed and no normal control was
available.

Several transcriptional states are strongly associated with individual
samples. Therefore, cell-level differential-expression and pathway results
are descriptive/exploratory and should not be interpreted as
patient-independent population-level effects.

Cells are not independent biological replicates.

CNV-like expression scores provide supportive malignancy evidence only and
are not equivalent to DNA-based CNV calling.
"""

    (out / "key_findings.md").write_text(findings)

    print("\n=== CELLFLOWX PHASE 12 COMPLETE ===")
    print(f"Samples summarized: {len(biological_summary)}")
    print(f"Pathway contrasts summarized: {len(pathway_summary)}")
    print("\nOutputs:")
    print("  biological_summary.csv")
    print("  state_pathway_summary.csv")
    print("  cellflowx_biological_synthesis.png")
    print("  key_findings.md")


if __name__ == "__main__":
    main()
