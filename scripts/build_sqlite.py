import sqlite3
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

DB_PATH = ROOT / "database" / "cellflowx.db"
SCHEMA_PATH = ROOT / "database" / "schema.sql"

MANIFEST = ROOT / "manifests" / "samplesheet.csv"
HET_DIR = ROOT / "results" / "heterogeneity"
CNV_DIR = ROOT / "results" / "cnv"
DE_DIR = ROOT / "results" / "differential_expression" / "phase10"
ENRICH_DIR = ROOT / "results" / "enrichment" / "phase11"


def load_manifest(conn):
    df = pd.read_csv(MANIFEST)

    out = pd.DataFrame({
        "sample_id": df["sample_id"],
        "geo_accession": df["geo_accession"],
        "condition": df["condition"],
        "tissue": df["tissue"],
        "organism": df["organism"],
        "assay": df["assay"],
        "n_cells": df["n_cells"],
        "n_genes": df["n_genes"],
    })

    out.to_sql(
        "samples",
        conn,
        if_exists="append",
        index=False
    )


def load_composition(conn):
    path = HET_DIR / "celltype_percentages_by_sample.csv"
    df = pd.read_csv(path)

    sample_col = (
        "sample"
        if "sample" in df.columns
        else "geo_accession"
    )

    long = df.melt(
        id_vars=sample_col,
        var_name="cell_state",
        value_name="percentage"
    )

    long = long.rename(
        columns={sample_col: "sample_id"}
    )

    long.to_sql(
        "celltype_composition",
        conn,
        if_exists="append",
        index=False
    )


def load_heterogeneity(conn):
    df = pd.read_csv(
        HET_DIR / "heterogeneity_summary.csv"
    )

    df = df.rename(
        columns={"sample": "sample_id"}
    )

    df.to_sql(
        "heterogeneity_summary",
        conn,
        if_exists="append",
        index=False
    )


def load_cnv(conn):
    df = pd.read_csv(
        CNV_DIR / "cnv_summary.csv"
    )

    rename_map = {
        "geo_accession": "sample_id",
        "cell_type": "cell_state",
    }

    df = df.rename(columns=rename_map)

    cols = [
        c for c in [
            "sample_id",
            "cell_state",
            "cnv_high_pct",
            "median_cnv_score"
        ]
        if c in df.columns
    ]

    df[cols].to_sql(
        "cnv_summary",
        conn,
        if_exists="append",
        index=False
    )


def load_de(conn):
    path = DE_DIR / "de_summary.csv"

    if not path.exists():
        candidates = list(DE_DIR.glob("*summary*.csv"))

        if not candidates:
            print("DE summary not found — skipped")
            return

        path = candidates[0]

    df = pd.read_csv(path)

    rename_map = {
        "state_A": "state_a",
        "state_B": "state_b",
        "cells_A": "cells_a",
        "cells_B": "cells_b",
        "up_in_A": "up_in_a",
        "up_in_B": "up_in_b",
    }

    df = df.rename(columns=rename_map)

    df.to_sql(
        "de_summary",
        conn,
        if_exists="append",
        index=False
    )


def load_enrichment(conn):
    df = pd.read_csv(
        ENRICH_DIR / "functional_enrichment_summary.csv"
    )

    df.to_sql(
        "enrichment_summary",
        conn,
        if_exists="append",
        index=False
    )


def main():
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)

    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)

    with open(SCHEMA_PATH) as handle:
        conn.executescript(handle.read())

    load_manifest(conn)
    load_composition(conn)
    load_heterogeneity(conn)
    load_cnv(conn)
    load_de(conn)
    load_enrichment(conn)

    conn.commit()

    print(f"SQLite database created: {DB_PATH}")

    for table in [
        "samples",
        "celltype_composition",
        "heterogeneity_summary",
        "cnv_summary",
        "de_summary",
        "enrichment_summary",
    ]:
        count = conn.execute(
            f"SELECT COUNT(*) FROM {table}"
        ).fetchone()[0]

        print(f"{table}: {count} rows")

    conn.close()


if __name__ == "__main__":
    main()
