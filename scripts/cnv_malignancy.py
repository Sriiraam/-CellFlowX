import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse
from pybiomart import Dataset


REFERENCE_TYPES = [
    "T/NK lymphocyte",
    "Macrophage",
    "Endothelial",
    "Fibroblast",
    "Activated fibroblast"
]

CANDIDATE_TYPES = [
    "Prostate epithelial - AR high",
    "Prostate epithelial - luminal",
    "Neuroendocrine-like",
    "Steroidogenic-like",
    "Epithelial-like - uncertain"
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--cell-scores", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--umap", required=True)
    p.add_argument("--heatmap", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    adata = sc.read_h5ad(args.input)

    if "cell_type" not in adata.obs:
        raise ValueError("cell_type annotation missing")

    if "gene_symbols" not in adata.var:
        raise ValueError("gene_symbols missing from adata.var")

    print(
        f"Input: {adata.n_obs:,} cells × "
        f"{adata.n_vars:,} genes"
    )

    # -------------------------------------------------
    # Obtain genomic positions from Ensembl
    # -------------------------------------------------

    print("\nRetrieving GRCh38 gene coordinates...")

    dataset = Dataset(
        name="hsapiens_gene_ensembl",
        host="http://www.ensembl.org"
    )

    genes = dataset.query(
        attributes=[
            "ensembl_gene_id",
            "chromosome_name",
            "start_position",
            "end_position"
        ]
    )

    genes.columns = [
        "ensembl_id",
        "chromosome",
        "start",
        "end"
    ]

    genes["ensembl_id"] = (
        genes["ensembl_id"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    genes = genes[
        genes["chromosome"].isin(
            [str(x) for x in range(1, 23)] + ["X", "Y"]
        )
    ].copy()

    genes = genes.drop_duplicates(
        "ensembl_id"
    )

    var_ids = (
        adata.var_names
        .astype(str)
        .str.split(".")
        .str[0]
    )

    coordinate_map = genes.set_index(
        "ensembl_id"
    )

    chromosome = []
    start = []

    for gene in var_ids:
        if gene in coordinate_map.index:
            row = coordinate_map.loc[gene]
            chromosome.append(str(row["chromosome"]))
            start.append(float(row["start"]))
        else:
            chromosome.append(np.nan)
            start.append(np.nan)

    adata.var["chromosome"] = chromosome
    adata.var["genomic_start"] = start

    usable = (
        adata.var["chromosome"].notna()
        & adata.var["genomic_start"].notna()
    )

    print(
        "Genes with genomic coordinates:",
        f"{usable.sum():,}"
    )

    if usable.sum() < 5000:
        raise ValueError(
            "Too few genes mapped to genomic coordinates."
        )

    # -------------------------------------------------
    # Reference baseline
    # -------------------------------------------------

    reference_mask = (
        adata.obs["cell_type"]
        .astype(str)
        .isin(REFERENCE_TYPES)
        .to_numpy()
    )

    candidate_mask = (
        adata.obs["cell_type"]
        .astype(str)
        .isin(CANDIDATE_TYPES)
        .to_numpy()
    )

    print(
        "Reference cells:",
        f"{reference_mask.sum():,}"
    )

    print(
        "Candidate epithelial/state cells:",
        f"{candidate_mask.sum():,}"
    )

    if reference_mask.sum() < 100:
        raise ValueError(
            "Insufficient non-malignant reference cells."
        )

    X = adata[:, usable].X

    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)

    X = X.tocsr()

    reference_mean = np.asarray(
        X[reference_mask].mean(axis=0)
    ).ravel()

    # -------------------------------------------------
    # Chromosome-level expression deviations
    # -------------------------------------------------

    var = adata.var.loc[usable].copy()

    chromosomes = [
        str(x) for x in range(1, 23)
    ] + ["X"]

    chromosome_scores = []

    for chrom in chromosomes:

        idx = np.where(
            var["chromosome"].astype(str).values == chrom
        )[0]

        if len(idx) < 20:
            continue

        chrom_expression = np.asarray(
            X[:, idx].mean(axis=1)
        ).ravel()

        chrom_reference = float(
            np.mean(reference_mean[idx])
        )

        deviation = (
            chrom_expression -
            chrom_reference
        )

        chromosome_scores.append(
            (chrom, deviation)
        )

    if not chromosome_scores:
        raise ValueError(
            "No chromosomes available for CNV scoring."
        )

    cnv_matrix = np.column_stack(
        [x[1] for x in chromosome_scores]
    )

    chrom_names = [
        x[0] for x in chromosome_scores
    ]

    # Mean absolute chromosome deviation
    cnv_score = np.mean(
        np.abs(cnv_matrix),
        axis=1
    )

    adata.obs["cnv_score"] = cnv_score

    # -------------------------------------------------
    # Reference-calibrated threshold
    # -------------------------------------------------

    ref_scores = cnv_score[
        reference_mask
    ]

    threshold = (
        np.median(ref_scores)
        + 3 * np.median(
            np.abs(
                ref_scores -
                np.median(ref_scores)
            )
        )
    )

    adata.obs["cnv_high"] = (
        adata.obs["cnv_score"] > threshold
    )

    adata.obs["malignancy_evidence"] = "not_assessed"

    adata.obs.loc[
        reference_mask,
        "malignancy_evidence"
    ] = "reference_like"

    adata.obs.loc[
        candidate_mask &
        adata.obs["cnv_high"].to_numpy(),
        "malignancy_evidence"
    ] = "CNV_high_candidate"

    adata.obs.loc[
        candidate_mask &
        ~adata.obs["cnv_high"].to_numpy(),
        "malignancy_evidence"
    ] = "CNV_low_candidate"

    # -------------------------------------------------
    # Tables
    # -------------------------------------------------

    cell_scores = adata.obs[
        [
            "geo_accession",
            "leiden",
            "cell_type",
            "cnv_score",
            "cnv_high",
            "malignancy_evidence"
        ]
    ].copy()

    cell_scores.to_csv(
        args.cell_scores
    )

    summary = (
        adata.obs
        .groupby(
            ["geo_accession", "cell_type"],
            observed=True
        )
        .agg(
            cells=("cnv_score", "size"),
            median_cnv_score=("cnv_score", "median"),
            mean_cnv_score=("cnv_score", "mean"),
            cnv_high_cells=("cnv_high", "sum")
        )
        .reset_index()
    )

    summary["cnv_high_pct"] = (
        summary["cnv_high_cells"]
        / summary["cells"]
        * 100
    ).round(2)

    summary.to_csv(
        args.summary,
        index=False
    )

    # -------------------------------------------------
    # UMAP
    # -------------------------------------------------

    sc.pl.umap(
        adata,
        color="cnv_score",
        show=False,
        title="CellFlowX — CNV Deviation Score"
    )

    plt.savefig(
        args.umap,
        dpi=220,
        bbox_inches="tight"
    )

    plt.close("all")

    # -------------------------------------------------
    # Chromosome heatmap
    # -------------------------------------------------

    candidate_indices = np.where(
        candidate_mask
    )[0]

    # Aggregate candidate CNV patterns by cell type.
    rows = []
    labels = []

    for cell_type in CANDIDATE_TYPES:

        mask = (
            adata.obs["cell_type"]
            .astype(str)
            .values == cell_type
        )

        if mask.sum() == 0:
            continue

        rows.append(
            np.median(
                cnv_matrix[mask],
                axis=0
            )
        )

        labels.append(cell_type)

    heat = np.vstack(rows)

    fig, ax = plt.subplots(
        figsize=(12, 5)
    )

    im = ax.imshow(
        heat,
        aspect="auto"
    )

    ax.set_xticks(
        range(len(chrom_names))
    )

    ax.set_xticklabels(
        chrom_names
    )

    ax.set_yticks(
        range(len(labels))
    )

    ax.set_yticklabels(
        labels
    )

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Candidate cell state")

    ax.set_title(
        "CellFlowX — Chromosome-Level Expression Deviations"
    )

    fig.colorbar(
        im,
        ax=ax,
        label="Expression deviation from reference"
    )

    fig.tight_layout()

    fig.savefig(
        args.heatmap,
        dpi=220,
        bbox_inches="tight"
    )

    plt.close(fig)

    # -------------------------------------------------
    # Metadata
    # -------------------------------------------------

    adata.uns["cnv_assessment"] = {
        "method":
            "chromosome-level expression deviation proxy",
        "reference_types":
            REFERENCE_TYPES,
        "candidate_types":
            CANDIDATE_TYPES,
        "threshold":
            float(threshold),
        "interpretation":
            "supportive malignancy evidence only; not DNA CNV calling"
    }

    adata.write_h5ad(
        args.output,
        compression="gzip"
    )

    print("\n=== CNV ASSESSMENT ===")
    print("Threshold:", round(float(threshold), 4))

    print("\nCNV-high cells:")
    print(
        adata.obs[
            "malignancy_evidence"
        ].value_counts()
    )

    print("\n=== PHASE 9 COMPLETE ===")

    for x in [
        args.output,
        args.cell_scores,
        args.summary,
        args.umap,
        args.heatmap
    ]:
        print(Path(x).resolve())


if __name__ == "__main__":
    main()
