import sqlite3
import pandas as pd
import streamlit as st

from utils.paths import PROJECT_ROOT

DB_PATH = PROJECT_ROOT / "database" / "cellflowx.db"


def get_connection():
    return sqlite3.connect(
        DB_PATH,
        check_same_thread=False
    )


@st.cache_data
def run_query(query, params=()):
    with get_connection() as conn:
        return pd.read_sql_query(
            query,
            conn,
            params=params
        )


@st.cache_data
def get_samples():
    return run_query(
        """
        SELECT *
        FROM samples
        ORDER BY geo_accession
        """
    )


@st.cache_data
def get_composition():
    return run_query(
        """
        SELECT
            sample_id,
            cell_state,
            percentage
        FROM celltype_composition
        ORDER BY sample_id, percentage DESC
        """
    )


@st.cache_data
def get_heterogeneity():
    return run_query(
        """
        SELECT *
        FROM heterogeneity_summary
        ORDER BY sample_id
        """
    )


@st.cache_data
def get_cnv_summary():
    return run_query(
        """
        SELECT *
        FROM cnv_summary
        ORDER BY sample_id, cnv_high_pct DESC
        """
    )


@st.cache_data
def get_de_summary():
    return run_query(
        """
        SELECT *
        FROM de_summary
        ORDER BY significant_genes DESC
        """
    )


@st.cache_data
def get_enrichment_summary():
    return run_query(
        """
        SELECT *
        FROM enrichment_summary
        ORDER BY top_adjusted_p
        """
    )
