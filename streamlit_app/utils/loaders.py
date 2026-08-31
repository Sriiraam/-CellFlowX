import pandas as pd
import streamlit as st


@st.cache_data
def load_csv(path):
    return pd.read_csv(path)


@st.cache_data
def load_table(path, index_col=None):
    return pd.read_csv(path, index_col=index_col)


@st.cache_data
def load_text(path):
    return path.read_text()


@st.cache_resource
def load_h5ad(path):
    import scanpy as sc
    return sc.read_h5ad(path)
