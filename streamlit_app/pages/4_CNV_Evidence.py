import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from utils.styles import apply_global_style, metric_card
from utils.paths import PROJECT_ROOT, CNV_DIR
from utils.loaders import load_csv, load_h5ad


st.set_page_config(
    page_title="CNV Evidence | CellFlowX",
    page_icon="🧬",
    layout="wide",
)

apply_global_style()

PALETTE = [
    "#6D4C91", "#9C6ADE", "#D77FA1", "#E89B6D",
    "#D4AD5C", "#6FA89A", "#75A6C9", "#A98CC1",
    "#C96F86", "#78A878", "#B7865B", "#8C7AA9",
]


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">CELLFLOWX • CNV-LIKE EVIDENCE</div>
<div class="hero-title">Expression-Based CNV Evidence</div>
<div class="hero-subtitle">
Explore chromosome-level expression deviation as supportive evidence
for tumor-associated cellular states across metastatic prostate cancer.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# LOAD
# --------------------------------------------------

CNV_CELLS_CSV = PROJECT_ROOT / "streamlit_app/data/cnv_cells.csv"
obs = load_csv(CNV_CELLS_CSV)

scores = load_csv(
    CNV_DIR / "cnv_cell_scores.csv"
)

summary = load_csv(
    CNV_DIR / "cnv_summary.csv"
)

# obs already loaded from deployment CSV
# UMAP coordinates already present in deployment CSV

def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


sample_col = find_col(
    scores,
    ["geo_accession", "sample", "sample_id"]
)

state_col = find_col(
    scores,
    ["cell_type", "cell_state", "annotation"]
)

score_col = find_col(
    scores,
    ["cnv_score", "score"]
)

label_col = find_col(
    scores,
    ["cnv_label", "label", "cnv_status"]
)


# --------------------------------------------------
# KPIs
# --------------------------------------------------

total_cells = len(scores)

if label_col:
    high_cells = (
        scores[label_col]
        .astype(str)
        .str.contains("CNV_high", case=False, na=False)
        .sum()
    )
else:
    high_cells = 0

candidate_states = (
    summary["cell_type"].nunique()
    if "cell_type" in summary.columns
    else len(summary)
)

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "CELLS SCORED",
        f"{total_cells:,}",
        "Chromosome-level expression"
    )

with c2:
    metric_card(
        "CNV-HIGH CELLS",
        f"{int(high_cells):,}",
        "Supportive candidates"
    )

with c3:
    metric_card(
        "STATES ASSESSED",
        str(candidate_states),
        "Tumor-associated populations"
    )

with c4:
    metric_card(
        "EVIDENCE TYPE",
        "Expression",
        "Not DNA-level CNV"
    )


# --------------------------------------------------
# SCORE DISTRIBUTION
# --------------------------------------------------

st.markdown("## CNV Score Landscape")

if score_col and state_col:

    selected_state = st.selectbox(
        "Cell state",
        ["All states"] + sorted(
            scores[state_col].astype(str).unique()
        )
    )

    plot_scores = scores.copy()

    if selected_state != "All states":
        plot_scores = plot_scores[
            plot_scores[state_col].astype(str)
            == selected_state
        ]

    fig_violin = px.violin(
        plot_scores,
        x=state_col,
        y=score_col,
        color=state_col,
        box=True,
        points=False,
        color_discrete_sequence=PALETTE,
        title="Distribution of CNV-like expression scores",
    )

    fig_violin.update_layout(
        height=560,
        xaxis_title="Cell state",
        yaxis_title="CNV-like score",
        showlegend=False,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        fig_violin,
        width="stretch"
    )


# --------------------------------------------------
# CNV SUMMARY HEATMAP
# --------------------------------------------------

st.markdown("## CNV-High Evidence by Tumor and Cell State")

if {
    "geo_accession",
    "cell_type",
    "cnv_high_pct"
}.issubset(summary.columns):

    heat = summary.pivot_table(
        index="cell_type",
        columns="geo_accession",
        values="cnv_high_pct",
        fill_value=0
    )

    fig_heat = px.imshow(
        heat,
        aspect="auto",
        text_auto=".1f",
        color_continuous_scale=[
            "#FBF8FD",
            "#E5D7EF",
            "#C1A4D4",
            "#8E63AD",
            "#5A376F",
        ],
        title="CNV-high candidate fraction",
    )

    fig_heat.update_layout(
        height=560,
        xaxis_title="Tumor",
        yaxis_title="Cell state",
        coloraxis_colorbar_title="CNV-high %",
        paper_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        fig_heat,
        width="stretch"
    )


# --------------------------------------------------
# UMAP EVIDENCE
# --------------------------------------------------

st.markdown("## Spatial View in Transcriptional Space")

obs_sample_col = find_col(
    obs,
    ["geo_accession", "sample", "sample_id"]
)

obs_state_col = find_col(
    obs,
    ["cell_type", "cell_state", "annotation"]
)

obs_score_col = find_col(
    obs,
    ["cnv_score", "score"]
)

if (
    "UMAP1" in obs.columns
    and "UMAP2" in obs.columns
    and obs_score_col
):

    fig_umap = px.scatter(
        obs,
        x="UMAP1",
        y="UMAP2",
        color=obs_score_col,
        hover_data=[
            c for c in [
                obs_sample_col,
                obs_state_col,
                obs_score_col
            ]
            if c is not None
        ],
        color_continuous_scale=[
            "#F7F2FA",
            "#D6BCE5",
            "#A879C5",
            "#6A4389",
        ],
        title="CNV-like score across the single-cell atlas",
        opacity=0.75,
    )

    fig_umap.update_traces(
        marker=dict(size=5)
    )

    fig_umap.update_layout(
        height=650,
        xaxis_title="UMAP 1",
        yaxis_title="UMAP 2",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        coloraxis_colorbar_title="CNV score",
    )

    fig_umap.update_xaxes(showgrid=False, zeroline=False)
    fig_umap.update_yaxes(showgrid=False, zeroline=False)

    st.plotly_chart(
        fig_umap,
        width="stretch"
    )


# --------------------------------------------------
# CANDIDATE EXPLORER
# --------------------------------------------------

st.markdown("## Candidate-State Explorer")

if {
    "geo_accession",
    "cell_type",
    "cnv_high_pct"
}.issubset(summary.columns):

    selected_candidate = st.selectbox(
        "Choose tumor-associated state",
        sorted(summary["cell_type"].astype(str).unique()),
        key="candidate_state"
    )

    candidate = summary[
        summary["cell_type"].astype(str)
        == selected_candidate
    ].copy()

    fig_bubble = px.scatter(
        candidate,
        x="geo_accession",
        y="cnv_high_pct",
        size="cnv_high_pct",
        color="geo_accession",
        color_discrete_sequence=[
            "#6D4C91",
            "#D77FA1",
            "#6FA89A"
        ],
        text=candidate[
            "cnv_high_pct"
        ].map(lambda x: f"{x:.1f}%"),
        size_max=75,
        title=f"{selected_candidate} • CNV-high evidence across tumors",
    )

    fig_bubble.update_traces(
        textposition="top center",
        marker=dict(
            line=dict(
                width=2,
                color="white"
            )
        )
    )

    fig_bubble.update_layout(
        height=430,
        xaxis_title="Tumor sample",
        yaxis_title="CNV-high cells (%)",
        showlegend=False,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        fig_bubble,
        width="stretch"
    )


# --------------------------------------------------
# STRONGEST EVIDENCE CARDS
# --------------------------------------------------

st.markdown("## Strongest CNV-Like Signals")

c1, c2, c3 = st.columns(3)

with c1:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">LUMINAL EPITHELIAL</div>
<div class="sample-state">Strongest signal</div>
<div class="sample-pct">40.1%</div>
<div class="metric-subtitle">
CNV-high candidates within GSM8848585 luminal cells.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c2:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">AR-HIGH EPITHELIAL</div>
<div class="sample-state">Substantial evidence</div>
<div class="sample-pct">24.2%</div>
<div class="metric-subtitle">
CNV-high candidates within GSM8848584 AR-high cells.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c3:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">OTHER STATES</div>
<div class="sample-state">Lower / mixed evidence</div>
<div class="metric-subtitle">
Steroidogenic-like and neuroendocrine-like states show lower
CNV-high fractions in this expression-based proxy.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


# --------------------------------------------------
# INTERPRETATION
# --------------------------------------------------

st.info(
    """
    **Interpretation:** chromosome-level expression deviation provides
    supportive evidence of CNV-like transcriptional alterations,
    particularly in AR-high and luminal prostate epithelial populations.
    This is an RNA-expression proxy and must not be interpreted as
    definitive DNA-level copy-number calling or proof of malignancy.
    """
)

with st.expander("Method summary"):

    st.markdown(
        """
        CellFlowX compares chromosome-level mean expression against
        reference-like stromal and immune populations. A per-cell
        expression-deviation score is calculated, and elevated scores
        are classified using a robust reference-based threshold.

        Reference populations include T/NK lymphocytes, macrophages,
        endothelial cells and fibroblast populations.

        Because transcription, technical noise and cell-state biology
        can all influence chromosome-level expression, the output is used
        strictly as **malignancy-supporting evidence**.
        """
    )


st.caption(
    "CellFlowX • CNV-like expression evidence • Supportive malignancy assessment"
)
