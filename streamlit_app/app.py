import streamlit as st
import plotly.express as px

from utils.styles import apply_global_style, metric_card
from utils.paths import HETEROGENEITY_DIR
from utils.loaders import load_csv
from utils.database import get_composition, get_heterogeneity


st.set_page_config(
    page_title="CellFlowX",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

apply_global_style()


# -----------------------------
# HERO
# -----------------------------

st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">SINGLE-CELL CANCER TRANSCRIPTOMICS</div>
<div class="hero-title">CellFlowX</div>
<div class="hero-subtitle">Interactive single-cell transcriptomic exploration of tumor heterogeneity in metastatic prostate cancer, powered by a reproducible Nextflow DSL2 workflow.</div>
</div>
""",
    unsafe_allow_html=True,
)


# -----------------------------
# LOAD DATA
# -----------------------------

summary = get_heterogeneity().rename(columns={"sample_id": "sample"})

composition_long = get_composition()

composition = composition_long.pivot(
    index="sample_id",
    columns="cell_state",
    values="percentage"
).reset_index().rename(columns={"sample_id": "sample"})

total_cells = int(summary["total_cells"].sum())
samples = len(summary)

state_columns = [
    c for c in composition.columns
    if c not in ["sample", "geo_accession"]
]

detected_states = len(state_columns)


# -----------------------------
# KPI CARDS
# -----------------------------

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "HIGH-QUALITY CELLS",
        f"{total_cells:,}",
        "Retained after QC"
    )

with c2:
    metric_card(
        "TUMOR SAMPLES",
        str(samples),
        "Metastatic tumors"
    )

with c3:
    metric_card(
        "CELL STATES",
        str(detected_states),
        "Annotated populations"
    )

with c4:
    metric_card(
        "PIPELINE",
        "DSL2",
        "Nextflow workflow"
    )


st.markdown("")


# -----------------------------
# OVERVIEW
# -----------------------------

st.markdown("## Study Overview")

left, right = st.columns([1.45, 1])

with left:
    st.markdown(
        """
        <div class="glass-card">

        <div class="card-title">
        Biological Objective
        </div>

        <div class="card-text">
        CellFlowX investigates how cellular composition and transcriptional
        states vary across metastatic prostate cancer tumors.

        The workflow integrates quality control, doublet detection,
        dimensionality reduction, Leiden clustering, marker-based annotation,
        tumor heterogeneity analysis, CNV-like expression evidence,
        differential expression and pathway enrichment.
        </div>

        </div>
        """,
        unsafe_allow_html=True,
    )


with right:
    st.markdown(
        """
        <div class="glass-card">

        <div class="card-title">
        Experimental Design
        </div>

        <div class="card-text">

        <b>Organism</b> — Homo sapiens<br>
        <b>Disease</b> — Metastatic prostate cancer<br>
        <b>Samples</b> — 3 metastatic tumors<br>
        <b>Assay</b> — 10x Genomics scRNA-seq<br>
        <b>Input</b> — Cell Ranger filtered matrices<br>
        <b>Control cohort</b> — None

        </div>

        </div>
        """,
        unsafe_allow_html=True,
    )


# -----------------------------
# INTERACTIVE COMPOSITION
# -----------------------------

st.markdown("## Tumor Landscape")

plot_df = composition.copy()

sample_col = (
    "sample"
    if "sample" in plot_df.columns
    else "geo_accession"
)

long_df = plot_df.melt(
    id_vars=sample_col,
    var_name="Cell state",
    value_name="Percentage"
)

fig = px.bar(
    long_df,
    x=sample_col,
    y="Percentage",
    color="Cell state",
    title="Interactive cellular composition across tumors",
)

fig.update_layout(
    height=470,
    margin=dict(l=10, r=10, t=55, b=10),
    legend_title_text="Cell state",
    xaxis_title="Tumor sample",
    yaxis_title="Cell composition (%)",
    hovermode="x unified",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(
    fig,
    width="stretch"
)


# -----------------------------
# DOMINANT BIOLOGY
# -----------------------------

st.markdown("## Dominant Tumor States")

cols = st.columns(3)

for col, (_, row) in zip(cols, summary.iterrows()):
    with col:
        html = f"""
<div class="sample-card">
<div class="sample-id">{row['sample']}</div>
<div class="sample-state">{row['dominant_state']}</div>
<div class="sample-pct">{row['dominant_state_pct']:.1f}%</div>
<div class="metric-subtitle">Dominant cellular state</div>
</div>
"""
        st.markdown(html, unsafe_allow_html=True)



# -----------------------------
# QUICK EXPLORER
# -----------------------------

st.markdown("## Quick Explorer")

selected_sample = st.selectbox(
    "Select tumor sample",
    summary["sample"].tolist(),
)

selected = summary[
    summary["sample"] == selected_sample
].iloc[0]

c1, c2, c3 = st.columns(3)

c1.metric(
    "Cells",
    f"{int(selected['total_cells']):,}"
)

c2.metric(
    "Detected states",
    int(selected["detected_states"])
)

c3.metric(
    "Dominant state",
    selected["dominant_state"]
)


# -----------------------------
# WORKFLOW
# -----------------------------

with st.expander("View CellFlowX analytical workflow"):

    st.markdown(
        """
        **10x matrices**
        → QC
        → Doublet detection
        → Filtering
        → Normalization
        → HVG selection
        → PCA
        → Neighbors
        → Leiden clustering
        → UMAP
        → Marker discovery
        → Cell-state annotation
        → Tumor heterogeneity
        → CNV-like expression evidence
        → Differential expression
        → GO enrichment
        → Biological synthesis
        """
    )


# -----------------------------
# INTERPRETATION NOTE
# -----------------------------

st.info(
    """
    **Interpretation note**

    CellFlowX analyzes three metastatic tumor samples without a normal-control
    cohort. Several transcriptional states are strongly associated with
    individual tumors, so differential-expression and pathway analyses are
    interpreted as exploratory cell-state comparisons rather than
    patient-level population inference.
    """
)


st.caption(
    "CellFlowX • Single-cell transcriptomics • Nextflow DSL2 • Scanpy • Streamlit"
)
