import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from utils.styles import apply_global_style, metric_card
from utils.paths import PROJECT_ROOT, QC_DIR
from utils.loaders import load_csv, load_h5ad


# --------------------------------------------------
# PAGE CONFIG
# --------------------------------------------------

st.set_page_config(
    page_title="Quality Control | CellFlowX",
    page_icon="🧬",
    layout="wide",
)

apply_global_style()


# --------------------------------------------------
# PATHS
# --------------------------------------------------

DEPLOY_DATA = PROJECT_ROOT / "streamlit_app/data"
QC_METRICS_CSV = DEPLOY_DATA / "qc_cell_metrics.csv"
QC_COUNTS_CSV = DEPLOY_DATA / "qc_counts.csv"


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">CELLFLOWX • QUALITY CONTROL</div>
<div class="hero-title">Sequencing Quality & Cell Filtering</div>
<div class="hero-subtitle">
Interactive assessment of library complexity, detected genes,
mitochondrial RNA content, sample-specific filtering and doublet removal.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

obs = load_csv(QC_METRICS_CSV)
qc_counts = load_csv(QC_COUNTS_CSV)

raw_cells = int(qc_counts.loc[0, "raw_cells"])
retained_cells = int(qc_counts.loc[0, "retained_cells"])
removed_cells = raw_cells - retained_cells
retention = retained_cells / raw_cells * 100


# --------------------------------------------------
# DETECT COLUMN NAMES
# --------------------------------------------------

def find_column(df, candidates):
    for col in candidates:
        if col in df.columns:
            return col
    return None


sample_col = find_column(
    obs,
    ["geo_accession", "sample", "sample_id"]
)

counts_col = find_column(
    obs,
    ["total_counts", "n_counts", "total_count"]
)

genes_col = find_column(
    obs,
    ["n_genes_by_counts", "n_genes", "genes_detected"]
)

mt_col = find_column(
    obs,
    ["pct_counts_mt", "pct_mt", "mitochondrial_percent"]
)


# --------------------------------------------------
# KPI CARDS
# --------------------------------------------------

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "RAW CELLS",
        f"{raw_cells:,}",
        "Before filtering"
    )

with c2:
    metric_card(
        "RETAINED CELLS",
        f"{retained_cells:,}",
        "High-quality cells"
    )

with c3:
    metric_card(
        "REMOVED CELLS",
        f"{removed_cells:,}",
        "QC + doublet exclusion"
    )

with c4:
    metric_card(
        "RETENTION",
        f"{retention:.1f}%",
        "Overall QC retention"
    )


st.markdown("")


# --------------------------------------------------
# QC OVERVIEW
# --------------------------------------------------

st.markdown("## QC Overview")

left, right = st.columns([1.35, 1])

with left:

    flow_df = pd.DataFrame(
        {
            "Stage": ["Raw cells", "Retained cells"],
            "Cells": [raw_cells, retained_cells],
        }
    )

    fig_flow = px.bar(
        flow_df,
        x="Stage",
        y="Cells",
        text="Cells",
        title="Cell retention through quality control",
    )

    fig_flow.update_traces(
        texttemplate="%{text:,}",
        textposition="outside",
    )

    fig_flow.update_layout(
        height=390,
        showlegend=False,
        margin=dict(l=10, r=10, t=55, b=10),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        yaxis_title="Number of cells",
        xaxis_title="",
    )

    st.plotly_chart(fig_flow, width="stretch")


with right:

    st.markdown(
        f"""
<div class="glass-card">
<div class="card-title">QC outcome</div>
<div class="card-text">
CellFlowX started with <b>{raw_cells:,}</b> cells across three metastatic
tumor samples.<br><br>

After sample-aware quality filtering and doublet exclusion,
<b>{retained_cells:,}</b> cells were retained.<br><br>

Overall retention was <b>{retention:.2f}%</b>, corresponding to
<b>{removed_cells:,}</b> excluded cells.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


# --------------------------------------------------
# SAMPLE SELECTOR
# --------------------------------------------------

st.markdown("## Interactive Cell-Level QC")

if sample_col:

    samples = ["All samples"] + sorted(
        obs[sample_col].astype(str).unique().tolist()
    )

    selected_sample = st.selectbox(
        "Tumor sample",
        samples,
    )

    if selected_sample == "All samples":
        plot_obs = obs.copy()
    else:
        plot_obs = obs[
            obs[sample_col].astype(str) == selected_sample
        ].copy()

else:
    selected_sample = "All samples"
    plot_obs = obs.copy()


# --------------------------------------------------
# COUNTS VS GENES
# --------------------------------------------------

if counts_col and genes_col:

    hover_data = []

    if mt_col:
        hover_data.append(mt_col)

    if sample_col:
        hover_data.append(sample_col)

    fig_scatter = px.scatter(
        plot_obs,
        x=counts_col,
        y=genes_col,
        color=sample_col if selected_sample == "All samples" else None,
        hover_data=hover_data,
        opacity=0.55,
        title="Library complexity: total counts vs detected genes",
    )

    fig_scatter.update_traces(
        marker=dict(size=5)
    )

    fig_scatter.update_layout(
        height=520,
        margin=dict(l=10, r=10, t=55, b=10),
        xaxis_title="Total counts",
        yaxis_title="Detected genes",
        legend_title="Tumor sample",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        fig_scatter,
        width="stretch"
    )

else:
    st.warning(
        "Counts/genes columns were not detected in the QC AnnData."
    )


# --------------------------------------------------
# DISTRIBUTIONS
# --------------------------------------------------

st.markdown("## QC Metric Distributions")

tab1, tab2, tab3 = st.tabs(
    [
        "Detected genes",
        "Total counts",
        "Mitochondrial RNA",
    ]
)


with tab1:

    if genes_col:

        fig = px.histogram(
            plot_obs,
            x=genes_col,
            nbins=70,
            color=sample_col
            if selected_sample == "All samples"
            else None,
            marginal="box",
            title="Genes detected per cell",
        )

        fig.update_layout(
            height=450,
            xaxis_title="Detected genes",
            yaxis_title="Cells",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )

        st.plotly_chart(fig, width="stretch")


with tab2:

    if counts_col:

        fig = px.histogram(
            plot_obs,
            x=counts_col,
            nbins=70,
            color=sample_col
            if selected_sample == "All samples"
            else None,
            marginal="box",
            title="Total RNA counts per cell",
        )

        fig.update_layout(
            height=450,
            xaxis_title="Total counts",
            yaxis_title="Cells",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )

        st.plotly_chart(fig, width="stretch")


with tab3:

    if mt_col:

        fig = px.histogram(
            plot_obs,
            x=mt_col,
            nbins=60,
            color=sample_col
            if selected_sample == "All samples"
            else None,
            marginal="box",
            title="Mitochondrial RNA percentage",
        )

        fig.update_layout(
            height=450,
            xaxis_title="Mitochondrial reads (%)",
            yaxis_title="Cells",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )

        st.plotly_chart(fig, width="stretch")


# --------------------------------------------------
# SAMPLE-WISE SUMMARY
# --------------------------------------------------

st.markdown("## Sample-Level Quality")

summary_path = QC_DIR / "qc_summary_by_sample.csv"

if summary_path.exists():

    summary = load_csv(summary_path)

    st.dataframe(
        summary,
        width="stretch",
        hide_index=True,
    )

else:
    st.info("Sample-level QC summary is unavailable.")


# --------------------------------------------------
# SAMPLE-SPECIFIC THRESHOLDS
# --------------------------------------------------

st.markdown("## Sample-Aware Filtering")

threshold_path = QC_DIR / "qc_threshold_evaluation.csv"

if threshold_path.exists():

    thresholds = load_csv(threshold_path)

    with st.expander(
        "View sample-specific QC threshold evaluation",
        expanded=False,
    ):
        st.dataframe(
            thresholds,
            width="stretch",
            hide_index=True,
        )


st.markdown(
    """
<div class="glass-card">
<div class="card-title">Why sample-specific thresholds?</div>
<div class="card-text">
The three tumors differ substantially in sequencing depth and cellular
composition. CellFlowX therefore avoids imposing one rigid global cutoff.
Gene-complexity and mitochondrial thresholds were evaluated independently
for each tumor before final filtering.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# DOUBLETS
# --------------------------------------------------

st.markdown("## Doublet Detection")

doublet_path = QC_DIR / "doublet_summary.csv"

if doublet_path.exists():

    doublets = load_csv(doublet_path)

    left, right = st.columns([1.3, 1])

    with left:

        st.dataframe(
            doublets,
            width="stretch",
            hide_index=True,
        )

    with right:

        st.markdown(
            """
<div class="glass-card">
<div class="card-title">Doublet strategy</div>
<div class="card-text">
Scrublet was executed independently for each tumor sample using an
expected doublet rate of 6%. Seven predicted doublet cells were
conservatively excluded from the final dataset.
<br><br>
Automatic global doublet-rate estimates were not used for biological
interpretation.
</div>
</div>
""",
            unsafe_allow_html=True,
        )


# --------------------------------------------------
# REMOVAL REASONS
# --------------------------------------------------

removal_path = QC_DIR / "qc_removal_reasons.csv"

if removal_path.exists():

    st.markdown("## Filtering Diagnostics")

    removal = load_csv(removal_path)

    with st.expander(
        "View QC removal diagnostics",
        expanded=False,
    ):
        st.dataframe(
            removal,
            width="stretch",
            hide_index=True,
        )


# --------------------------------------------------
# QC INTERPRETATION
# --------------------------------------------------

st.markdown("## QC Interpretation")

c1, c2, c3 = st.columns(3)

with c1:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">LIBRARY COMPLEXITY</div>
<div class="sample-state">Strong retained signal</div>
<div class="metric-subtitle">
Most cells retain substantial transcript and gene complexity after filtering.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


with c2:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">MITOCHONDRIAL QC</div>
<div class="sample-state">Sample-aware</div>
<div class="metric-subtitle">
Mitochondrial content varies between tumors and was handled using
sample-specific thresholds.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


with c3:

    st.markdown(
        f"""
<div class="sample-card">
<div class="sample-id">FINAL DATASET</div>
<div class="sample-state">{retained_cells:,} cells</div>
<div class="metric-subtitle">
High-quality sparse expression matrix carried forward into downstream analysis.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


# --------------------------------------------------
# IMPORTANT DATASET NOTE
# --------------------------------------------------

st.info(
    """
    **Dataset-specific QC note:** canonical ribosomal RPL/RPS genes are absent
    from the supplied processed feature matrices. Ribosomal percentage was
    therefore not fabricated or used as a QC filtering criterion.
    """
)


st.caption(
    "CellFlowX • Quality Control • Sample-aware filtering • Scrublet"
)
