import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from utils.styles import apply_global_style, metric_card
from utils.paths import PROJECT_ROOT, ANNOTATION_DIR
from utils.loaders import load_csv, load_h5ad


st.set_page_config(
    page_title="Cell Atlas | CellFlowX",
    page_icon="🧬",
    layout="wide",
)

apply_global_style()

PALETTE = [
    "#6D4C91", "#9C6ADE", "#D77FA1", "#E89B6D",
    "#D4AD5C", "#6FA89A", "#75A6C9", "#A98CC1",
    "#C96F86", "#78A878", "#B7865B", "#8C7AA9",
]


st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">CELLFLOWX • CELL ATLAS</div>
<div class="hero-title">Interactive Single-Cell Atlas</div>
<div class="hero-subtitle">
Explore tumor-specific transcriptional neighborhoods, annotated cell states,
cluster structure and marker programs across metastatic prostate cancer.
</div>
</div>
""",
    unsafe_allow_html=True,
)


ATLAS_CSV = PROJECT_ROOT / "streamlit_app/data/cell_atlas.csv"
obs = load_csv(ATLAS_CSV)

# UMAP coordinates already included in deployment CSV
plot_df = obs.copy()


def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


sample_col = find_col(
    plot_df,
    ["geo_accession", "sample", "sample_id"]
)

cluster_col = find_col(
    plot_df,
    ["leiden", "cluster", "clusters"]
)

celltype_col = find_col(
    plot_df,
    [
        "cell_type",
        "celltype",
        "cell_state",
        "annotation",
        "cell_type_annotation",
    ]
)

if celltype_col is None:
    st.error("Cell-state annotation column was not detected.")
    st.stop()


# KPIs
c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card("CELLS", f"{len(plot_df):,}", "QC-passed cells")

with c2:
    metric_card(
        "CELL STATES",
        plot_df[celltype_col].nunique(),
        "Annotated populations"
    )

with c3:
    metric_card(
        "LEIDEN CLUSTERS",
        plot_df[cluster_col].nunique() if cluster_col else "-",
        "Transcriptomic neighborhoods"
    )

with c4:
    metric_card(
        "TUMORS",
        plot_df[sample_col].nunique() if sample_col else "-",
        "Metastatic samples"
    )


st.markdown("## Atlas Explorer")

f1, f2, f3 = st.columns([1, 1.4, 1])

with f1:
    sample_options = ["All tumors"] + sorted(
        plot_df[sample_col].astype(str).unique()
    )

    selected_sample = st.selectbox(
        "Tumor sample",
        sample_options
    )

with f2:
    state_options = ["All cell states"] + sorted(
        plot_df[celltype_col].astype(str).unique()
    )

    selected_state = st.selectbox(
        "Cell state",
        state_options
    )

with f3:
    color_mode = st.selectbox(
        "Color UMAP by",
        ["Cell state", "Tumor sample", "Leiden cluster"]
    )


filtered = plot_df.copy()

if selected_sample != "All tumors":
    filtered = filtered[
        filtered[sample_col].astype(str) == selected_sample
    ]

if selected_state != "All cell states":
    filtered = filtered[
        filtered[celltype_col].astype(str) == selected_state
    ]


if color_mode == "Cell state":
    color_col = celltype_col
elif color_mode == "Tumor sample":
    color_col = sample_col
else:
    color_col = cluster_col


# --------------------------------------------------
# UMAP
# --------------------------------------------------

fig_umap = px.scatter(
    filtered,
    x="UMAP1",
    y="UMAP2",
    color=color_col,
    color_discrete_sequence=PALETTE,
    hover_data=[sample_col, celltype_col, cluster_col],
    opacity=0.78,
    title=f"Interactive UMAP • {len(filtered):,} cells",
)

fig_umap.update_traces(
    marker=dict(size=5)
)

fig_umap.update_layout(
    height=680,
    margin=dict(l=10, r=10, t=60, b=10),
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

fig_umap.update_xaxes(showgrid=False, zeroline=False)
fig_umap.update_yaxes(showgrid=False, zeroline=False)

st.plotly_chart(fig_umap, width="stretch")


# --------------------------------------------------
# DONUT + HEATMAP
# --------------------------------------------------

st.markdown("## Cellular Architecture")

left, right = st.columns([1, 1.45])


with left:

    counts = (
        plot_df[celltype_col]
        .value_counts()
        .reset_index()
    )

    counts.columns = ["Cell state", "Cells"]

    fig_donut = px.pie(
        counts,
        names="Cell state",
        values="Cells",
        hole=0.58,
        color_discrete_sequence=PALETTE,
        title="Cell-state composition"
    )

    fig_donut.update_traces(
        textposition="inside",
        textinfo="percent"
    )

    fig_donut.update_layout(
        height=520,
        margin=dict(l=10, r=10, t=55, b=10),
        legend_title="Cell state",
        paper_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(fig_donut, width="stretch")


with right:

    cluster_state = pd.crosstab(
        plot_df[cluster_col].astype(str),
        plot_df[celltype_col].astype(str),
        normalize="index"
    ) * 100

    fig_heat = px.imshow(
        cluster_state,
        aspect="auto",
        text_auto=".0f",
        color_continuous_scale=[
            "#F8F3FC",
            "#D8C2EA",
            "#9B78BE",
            "#5B3A78",
        ],
        title="Cluster × cell-state identity map"
    )

    fig_heat.update_layout(
        height=520,
        xaxis_title="Cell state",
        yaxis_title="Leiden cluster",
        coloraxis_colorbar_title="%",
        paper_bgcolor="rgba(0,0,0,0)",
    )

    st.plotly_chart(fig_heat, width="stretch")


# --------------------------------------------------
# LOLLIPOP DISTRIBUTION
# --------------------------------------------------

st.markdown("## Selected Population Across Tumors")

selected_population = st.selectbox(
    "Choose a population",
    sorted(plot_df[celltype_col].astype(str).unique()),
    key="lollipop_state"
)

state_df = plot_df[
    plot_df[celltype_col].astype(str) == selected_population
]

distribution = (
    state_df[sample_col]
    .value_counts()
    .rename_axis("Tumor")
    .reset_index(name="Cells")
)

all_samples = pd.DataFrame(
    {"Tumor": sorted(plot_df[sample_col].astype(str).unique())}
)

distribution = all_samples.merge(
    distribution,
    on="Tumor",
    how="left"
).fillna({"Cells": 0})


fig_lollipop = go.Figure()

for _, row in distribution.iterrows():

    fig_lollipop.add_trace(
        go.Scatter(
            x=[row["Tumor"], row["Tumor"]],
            y=[0, row["Cells"]],
            mode="lines",
            line=dict(width=4),
            showlegend=False,
            hoverinfo="skip"
        )
    )

fig_lollipop.add_trace(
    go.Scatter(
        x=distribution["Tumor"],
        y=distribution["Cells"],
        mode="markers+text",
        marker=dict(
            size=22,
            color=["#6D4C91", "#D77FA1", "#6FA89A"][:len(distribution)]
        ),
        text=distribution["Cells"].astype(int),
        textposition="top center",
        hovertemplate="<b>%{x}</b><br>Cells: %{y}<extra></extra>",
        showlegend=False
    )
)

fig_lollipop.update_layout(
    title=f"{selected_population} distribution across tumors",
    height=430,
    xaxis_title="Tumor sample",
    yaxis_title="Cells",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(fig_lollipop, width="stretch")


# --------------------------------------------------
# MARKER BUBBLE EXPLORER
# --------------------------------------------------

st.markdown("## Marker Bubble Explorer")

marker_path = ANNOTATION_DIR / "top_cluster_markers_symbols.csv"

if not marker_path.exists():
    marker_path = ANNOTATION_DIR / "top_cluster_markers.csv"


if marker_path.exists():

    markers = load_csv(marker_path)

    marker_cluster_col = find_col(
        markers,
        ["cluster", "group", "leiden"]
    )

    marker_gene_col = find_col(
        markers,
        ["gene_symbol", "gene", "names", "symbol"]
    )

    marker_score_col = find_col(
        markers,
        ["score", "scores", "logfoldchanges", "logFC"]
    )

    if (
        marker_cluster_col
        and marker_gene_col
        and marker_score_col
    ):

        bubble = markers.copy()

        bubble = bubble.groupby(
            marker_cluster_col,
            group_keys=False
        ).head(5)

        bubble["bubble_size"] = (
            bubble[marker_score_col].abs()
        )

        fig_bubble = px.scatter(
            bubble,
            x=marker_cluster_col,
            y=marker_gene_col,
            size="bubble_size",
            color=marker_cluster_col,
            color_discrete_sequence=PALETTE,
            hover_data=[marker_score_col],
            title="Top marker programs across Leiden clusters"
        )

        fig_bubble.update_layout(
            height=620,
            xaxis_title="Leiden cluster",
            yaxis_title="Marker gene",
            showlegend=False,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )

        st.plotly_chart(
            fig_bubble,
            width="stretch"
        )

    else:

        st.dataframe(
            markers,
            width="stretch",
            hide_index=True
        )


# --------------------------------------------------
# BIOLOGICAL INSIGHTS
# --------------------------------------------------

st.markdown("## Biological Insights")

c1, c2, c3 = st.columns(3)

with c1:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">EPITHELIAL DIVERSITY</div>
<div class="sample-state">AR-high & luminal states</div>
<div class="metric-subtitle">
Prostate epithelial populations occupy distinct transcriptional regions,
supporting substantial epithelial-state heterogeneity.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c2:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">NEUROENDOCRINE PROGRAMS</div>
<div class="sample-state">Distinct transcriptional state</div>
<div class="metric-subtitle">
Neuroendocrine-like clusters show strong separation and characteristic
ASCL1/INSM1-associated expression programs.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c3:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">MICROENVIRONMENT</div>
<div class="sample-state">Stromal & immune diversity</div>
<div class="metric-subtitle">
Fibroblast, macrophage, endothelial and lymphocyte populations vary
substantially across tumors.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


st.info(
    """
    UMAP represents transcriptional similarity, not physical tumor location,
    developmental trajectory or proof of malignant lineage.
    """
)

st.caption(
    "CellFlowX • Interactive Cell Atlas • UMAP • Cell states • Marker programs"
)