import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from utils.styles import apply_global_style, metric_card
from utils.paths import HETEROGENEITY_DIR
from utils.loaders import load_csv
from utils.database import get_composition, get_heterogeneity


st.set_page_config(
    page_title="Tumor Heterogeneity | CellFlowX",
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
<div class="hero-eyebrow">CELLFLOWX • INTER-TUMOR HETEROGENEITY</div>
<div class="hero-title">Three Tumors. Three Cellular Ecosystems.</div>
<div class="hero-subtitle">
Compare cellular composition across metastatic prostate cancer samples
and explore the tumor-specific states driving inter-tumor heterogeneity.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# LOAD
# --------------------------------------------------

summary = get_heterogeneity().rename(columns={"sample_id": "sample"})

composition_long = get_composition()

composition = composition_long.pivot(
    index="sample_id",
    columns="cell_state",
    values="percentage"
).reset_index().rename(columns={"sample_id": "sample"})

sample_col = (
    "sample"
    if "sample" in composition.columns
    else "geo_accession"
)

state_cols = [
    c for c in composition.columns
    if c != sample_col
]


# --------------------------------------------------
# KPIs
# --------------------------------------------------

total_cells = int(summary["total_cells"].sum())

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "TUMORS",
        len(summary),
        "Independent metastatic samples"
    )

with c2:
    metric_card(
        "CELLS PROFILED",
        f"{total_cells:,}",
        "After QC"
    )

with c3:
    metric_card(
        "CELL STATES",
        len(state_cols),
        "Atlas-wide populations"
    )

with c4:
    metric_card(
        "MAX DOMINANCE",
        f"{summary['dominant_state_pct'].max():.1f}%",
        "Largest tumor-specific state"
    )


# --------------------------------------------------
# 100% COMPOSITION
# --------------------------------------------------

st.markdown("## Cellular Composition Fingerprint")

long = composition.melt(
    id_vars=sample_col,
    var_name="Cell state",
    value_name="Percentage"
)

fig = px.bar(
    long,
    x=sample_col,
    y="Percentage",
    color="Cell state",
    color_discrete_sequence=PALETTE,
    title="Cell-state composition of each tumor",
)

fig.update_layout(
    barmode="stack",
    height=560,
    yaxis_range=[0, 100],
    xaxis_title="Tumor sample",
    yaxis_title="Cell composition (%)",
    legend_title="Cell state",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

fig.update_traces(
    hovertemplate=(
        "<b>%{x}</b><br>"
        "%{fullData.name}<br>"
        "%{y:.2f}%"
        "<extra></extra>"
    )
)

st.plotly_chart(fig, width="stretch")


# --------------------------------------------------
# HEATMAP + RADAR
# --------------------------------------------------

st.markdown("## Heterogeneity Explorer")

left, right = st.columns([1.35, 1])


# HEATMAP
with left:

    heat = composition.set_index(sample_col).T

    fig_heat = px.imshow(
        heat,
        aspect="auto",
        text_auto=".1f",
        color_continuous_scale=[
            "#FBF8FD",
            "#E2D3EE",
            "#B99AD0",
            "#7A5599",
            "#4B2E63",
        ],
        title="Cell-state enrichment map",
    )

    fig_heat.update_layout(
        height=610,
        xaxis_title="Tumor",
        yaxis_title="Cell state",
        coloraxis_colorbar_title="%",
        paper_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(
        fig_heat,
        width="stretch"
    )


# RADAR
with right:

    selected_tumor = st.selectbox(
        "Explore tumor profile",
        composition[sample_col].tolist()
    )

    row = composition[
        composition[sample_col] == selected_tumor
    ].iloc[0]

    values = [
        float(row[state])
        for state in state_cols
    ]

    fig_radar = go.Figure()

    fig_radar.add_trace(
        go.Scatterpolar(
            r=values + [values[0]],
            theta=state_cols + [state_cols[0]],
            fill="toself",
            line=dict(
                color="#6D4C91",
                width=3
            ),
            fillcolor="rgba(156,106,222,0.22)",
            hovertemplate=(
                "<b>%{theta}</b><br>"
                "%{r:.2f}%"
                "<extra></extra>"
            ),
        )
    )

    fig_radar.update_layout(
        title=f"{selected_tumor} cellular fingerprint",
        height=610,
        showlegend=False,
        paper_bgcolor="rgba(0,0,0,0)",
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            radialaxis=dict(
                visible=True,
                gridcolor="#DED3E7"
            )
        ),
        margin=dict(
            l=70,
            r=70,
            t=80,
            b=70
        )
    )

    st.plotly_chart(
        fig_radar,
        width="stretch"
    )


# --------------------------------------------------
# DOMINANT STATE CARDS
# --------------------------------------------------

st.markdown("## Dominant Biology by Tumor")

cols = st.columns(3)

for col, (_, row) in zip(
    cols,
    summary.iterrows()
):
    with col:

        st.markdown(
            f"""
<div class="sample-card">
<div class="sample-id">{row['sample']}</div>
<div class="sample-state">{row['dominant_state']}</div>
<div class="sample-pct">{row['dominant_state_pct']:.1f}%</div>
<div class="metric-subtitle">
dominant cellular population
</div>
<br>
<div class="card-text">
{int(row['total_cells']):,} cells •
{int(row['detected_states'])} detected states
</div>
</div>
""",
            unsafe_allow_html=True,
        )


# --------------------------------------------------
# STATE-CENTRIC EXPLORER
# --------------------------------------------------

st.markdown("## Which Tumor Drives Each Cell State?")

selected_state = st.selectbox(
    "Select cellular state",
    state_cols,
    key="state_driver"
)

state_distribution = composition[
    [sample_col, selected_state]
].copy()

state_distribution.columns = [
    "Tumor",
    "Percentage"
]


# Use bubble representation rather than another bar
fig_bubble = px.scatter(
    state_distribution,
    x="Tumor",
    y="Percentage",
    size="Percentage",
    color="Tumor",
    color_discrete_sequence=[
        "#6D4C91",
        "#D77FA1",
        "#6FA89A"
    ],
    text=state_distribution[
        "Percentage"
    ].map(lambda x: f"{x:.1f}%"),
    title=f"{selected_state} across tumors",
    size_max=70,
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
    yaxis_title="Cell composition (%)",
    showlegend=False,
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(
    fig_bubble,
    width="stretch"
)


# --------------------------------------------------
# TOP 3 STATES
# --------------------------------------------------

st.markdown("## Tumor Signatures")

tumor_cols = st.columns(3)

for col, (_, row) in zip(
    tumor_cols,
    composition.iterrows()
):

    tumor = row[sample_col]

    top = (
        row[state_cols]
        .astype(float)
        .sort_values(ascending=False)
        .head(3)
    )

    with col:

        st.markdown(
            f"### {tumor}"
        )

        for rank, (state, pct) in enumerate(
            top.items(),
            start=1
        ):

            st.progress(
                min(float(pct) / 100, 1.0),
                text=f"{rank}. {state} — {pct:.1f}%"
            )


# --------------------------------------------------
# INTERPRETATION
# --------------------------------------------------

st.markdown("## Biological Interpretation")

c1, c2, c3 = st.columns(3)

with c1:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">GSM8848584</div>
<div class="sample-state">AR-high / steroidogenic</div>
<div class="metric-subtitle">
The tumor is dominated by steroidogenic-like and AR-high
prostate epithelial transcriptional programs.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


with c2:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">GSM8848585</div>
<div class="sample-state">Neuroendocrine-rich</div>
<div class="metric-subtitle">
A large neuroendocrine-like compartment accompanies luminal
and proliferative cellular states.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


with c3:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">GSM8848586</div>
<div class="sample-state">Stromal / myeloid-rich</div>
<div class="metric-subtitle">
Activated fibroblasts dominate this sample together with
a substantial macrophage compartment.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


st.info(
    """
    **Main finding:** the three metastatic tumors occupy markedly different
    cellular ecosystems. This supports strong inter-tumor heterogeneity,
    while the small three-sample cohort means these patterns should remain
    descriptive rather than population-level estimates.
    """
)

st.caption(
    "CellFlowX • Inter-Tumor Heterogeneity • Composition • Cellular ecosystems"
)
