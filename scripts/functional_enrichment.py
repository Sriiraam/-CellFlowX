import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--de-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def find_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


def safe_name(x):
    return (
        str(x)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("-", "_")
    )


def run_enrichment(genes, label, output_dir):
    genes = [
        str(g).strip()
        for g in genes
        if pd.notna(g) and str(g).strip() not in {"", "nan", "None"}
    ]
    genes = list(dict.fromkeys(genes))

    print(f"\n{label}: {len(genes):,} genes")

    if len(genes) < 10:
        print("Skipping: fewer than 10 genes.")
        return None

    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets="GO_Biological_Process_2023",
            organism="human",
            outdir=None,
            cutoff=0.05,
        )
    except Exception as e:
        print(f"Enrichment failed for {label}: {e}")
        return None

    if enr.results is None or enr.results.empty:
        print("No significant pathways.")
        return None

    result = enr.results.copy()

    if "Adjusted P-value" in result.columns:
        result = result[result["Adjusted P-value"] < 0.05].copy()

    if result.empty:
        print("No pathways with adjusted P < 0.05.")
        return None

    result = result.sort_values("Adjusted P-value")

    stem = safe_name(label)

    csv_path = output_dir / f"{stem}_GO_BP.csv"
    result.to_csv(csv_path, index=False)

    top = result.head(15).copy()

    top["minus_log10_fdr"] = -np.log10(
        top["Adjusted P-value"].clip(lower=1e-300)
    )

    # Reverse so most significant appears at top.
    top = top.iloc[::-1]

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.barh(
        top["Term"],
        top["minus_log10_fdr"]
    )

    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_ylabel("")
    ax.set_title(label.replace("_", " ") + " — GO Biological Process")

    plt.tight_layout()

    plot_path = output_dir / f"{stem}_GO_BP.png"
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Significant pathways: {len(result):,}")
    print(f"Saved: {csv_path}")
    print(f"Saved: {plot_path}")

    return {
        "gene_set": label,
        "input_genes": len(genes),
        "significant_pathways": len(result),
        "top_pathway": result.iloc[0]["Term"],
        "top_adjusted_p": result.iloc[0]["Adjusted P-value"],
    }


def main():
    args = parse_args()

    de_dir = Path(args.de_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(de_dir.glob("*_DE.csv"))

    if not files:
        # Fall back to all CSVs except summaries.
        files = [
            p for p in sorted(de_dir.glob("*.csv"))
            if "summary" not in p.name.lower()
            and "representation" not in p.name.lower()
        ]

    if not files:
        raise FileNotFoundError(
            f"No Phase 10 DE result CSV files found in {de_dir}"
        )

    print("=== CELLFLOWX PHASE 11 ===")
    print(f"DE directory: {de_dir}")
    print(f"DE files found: {len(files)}")

    summaries = []

    for file in files:
        print(f"\nReading: {file.name}")

        df = pd.read_csv(file)

        gene_col = find_column(
            df,
            [
                "gene_symbol",
                "gene_symbols",
                "symbol",
                "Gene",
                "gene",
                "names",
            ],
        )

        lfc_col = find_column(
            df,
            [
                "logfoldchanges",
                "log2FoldChange",
                "logFC",
                "log2fc",
            ],
        )

        padj_col = find_column(
            df,
            [
                "pvals_adj",
                "padj",
                "adjusted_pvalue",
                "Adjusted P-value",
            ],
        )

        if gene_col is None or lfc_col is None or padj_col is None:
            print(
                "Skipping because required columns were not found.\n"
                f"Columns: {list(df.columns)}"
            )
            continue

        # If names are Ensembl IDs but a symbol column exists, prefer symbol.
        if gene_col == "names":
            symbol_col = find_column(
                df,
                ["gene_symbol", "gene_symbols", "symbol"]
            )
            if symbol_col is not None:
                gene_col = symbol_col

        sig = df[
            (pd.to_numeric(df[padj_col], errors="coerce") < 0.05)
            &
            (pd.to_numeric(df[lfc_col], errors="coerce").abs() >= 1)
        ].copy()

        if sig.empty:
            print("No significant genes after |logFC| >= 1 and FDR < 0.05.")
            continue

        comparison = file.stem
        comparison = comparison.replace("_DE", "")

        up_a = sig[
            pd.to_numeric(sig[lfc_col], errors="coerce") >= 1
        ][gene_col]

        up_b = sig[
            pd.to_numeric(sig[lfc_col], errors="coerce") <= -1
        ][gene_col]

        a_result = run_enrichment(
            up_a,
            f"{comparison}_up_A",
            output_dir,
        )

        b_result = run_enrichment(
            up_b,
            f"{comparison}_up_B",
            output_dir,
        )

        if a_result:
            summaries.append(a_result)

        if b_result:
            summaries.append(b_result)

    if not summaries:
        raise RuntimeError(
            "Phase 11 finished without producing any enrichment results. "
            "Check the DE filenames/columns and gene symbols."
        )

    summary = pd.DataFrame(summaries)

    summary_file = output_dir / "functional_enrichment_summary.csv"
    summary.to_csv(summary_file, index=False)

    print("\n=== PHASE 11 COMPLETE ===")
    print(summary.to_string(index=False))
    print(f"\nSaved summary: {summary_file}")


if __name__ == "__main__":
    main()
