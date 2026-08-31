import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from utils.styles import apply_global_style, metric_card
from utils.paths import DE_DIR
from utils.loaders import load_csv


st.set_page_config(
    page_title="Differential Expression | CellFlowX",
    page_icon="🧬",
    layout="wide",
)

apply_global_style()

COLORS = {
    "Up in A": "#6D4C91",
    "Up in B": "#D77FA1",
    "Not significant": "#CFC5D8",
}


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">CELLFLOWX • DIFFERENTIAL EXPRESSION</div>
<div class="hero-title">Tumor-State Gene Programs</div>
<div class="hero-subtitle">
Interactively explore transcriptional differences between major
tumor-associated cellular states and identify genes defining each program.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# FIND DE FILES
# --------------------------------------------------

de_files = sorted(DE_DIR.glob("*_DE.csv"))

if not de_files:
    st.error("No Phase 10 differential-expression files were found.")
    st.stop()


comparisons = {
    f.stem.replace("_DE", ""): f
    for f in de_files
}


# --------------------------------------------------
# COMPARISON SELECTOR
# --------------------------------------------------

st.markdown("## Comparison Explorer")

selected_comparison = st.selectbox(
    "Choose tumor-state comparison",
    list(comparisons.keys()),
    format_func=lambda x: x.replace("_", " "),
)

df = load_csv(comparisons[selected_comparison])


def find_col(data, candidates):
    for c in candidates:
        if c in data.columns:
            return c
    return None


gene_col = find_col(
    df,
    ["gene", "gene_symbol", "names", "symbol"]
)

lfc_col = find_col(
    df,
    ["logfoldchanges", "logFC", "log2FoldChange"]
)

padj_col = find_col(
    df,
    ["pvals_adj", "padj", "adjusted_pvalue"]
)


if not all([gene_col, lfc_col, padj_col]):
    st.error(
        "Could not identify gene, log-fold-change and adjusted-p columns."
    )
    st.stop()


# --------------------------------------------------
# DEFINE SIGNIFICANCE
# --------------------------------------------------

df = df.copy()

df[lfc_col] = pd.to_numeric(
    df[lfc_col],
    errors="coerce"
)

df[padj_col] = pd.to_numeric(
    df[padj_col],
    errors="coerce"
)

df = df.dropna(
    subset=[lfc_col, padj_col]
)

df["neglog10_padj"] = -np.log10(
    df[padj_col].clip(lower=1e-300)
)

df["Direction"] = "Not significant"

df.loc[
    (df[padj_col] < 0.05) &
    (df[lfc_col] >= 1),
    "Direction"
] = "Up in A"

df.loc[
    (df[padj_col] < 0.05) &
    (df[lfc_col] <= -1),
    "Direction"
] = "Up in B"


sig = df[
    df["Direction"] != "Not significant"
]

up_a = int(
    (df["Direction"] == "Up in A").sum()
)

up_b = int(
    (df["Direction"] == "Up in B").sum()
)


# --------------------------------------------------
# KPIs
# --------------------------------------------------

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "GENES TESTED",
        f"{len(df):,}",
        "Expression features"
    )

with c2:
    metric_card(
        "SIGNIFICANT",
        f"{len(sig):,}",
        "FDR < 0.05 & |logFC| ≥ 1"
    )

with c3:
    metric_card(
        "UP IN STATE A",
        f"{up_a:,}",
        "Positive logFC"
    )

with c4:
    metric_card(
        "UP IN STATE B",
        f"{up_b:,}",
        "Negative logFC"
    )


# --------------------------------------------------
# INTERACTIVE VOLCANO
# --------------------------------------------------

st.markdown("## Interactive Volcano")

fig = px.scatter(
    df,
    x=lfc_col,
    y="neglog10_padj",
    color="Direction",
    color_discrete_map=COLORS,
    hover_name=gene_col,
    hover_data={
        lfc_col: ":.2f",
        padj_col: ":.2e",
        "neglog10_padj": False,
    },
    opacity=0.65,
    title=selected_comparison.replace("_", " "),
)

fig.add_vline(
    x=-1,
    line_dash="dash",
    line_color="#9B8AA8"
)

fig.add_vline(
    x=1,
    line_dash="dash",
    line_color="#9B8AA8"
)

fig.add_hline(
    y=-np.log10(0.05),
    line_dash="dash",
    line_color="#9B8AA8"
)

fig.update_traces(
    marker=dict(size=6)
)

fig.update_layout(
    height=650,
    xaxis_title="Log fold change",
    yaxis_title="−log10 adjusted P-value",
    legend_title="DE direction",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(fig, width="stretch")


# --------------------------------------------------
# GENE SEARCH
# --------------------------------------------------

st.markdown("## Gene Inspector")

genes = sorted(
    df[gene_col]
    .dropna()
    .astype(str)
    .unique()
)

selected_gene = st.selectbox(
    "Search gene",
    genes
)

gene_row = df[
    df[gene_col].astype(str)
    == selected_gene
].iloc[0]


g1, g2, g3 = st.columns(3)

g1.metric(
    "Log fold change",
    f"{gene_row[lfc_col]:.3f}"
)

g2.metric(
    "Adjusted P-value",
    f"{gene_row[padj_col]:.2e}"
)

g3.metric(
    "Direction",
    gene_row["Direction"]
)


# --------------------------------------------------
# DIRECTIONAL TOP GENES
# --------------------------------------------------

st.markdown("## Leading Differential Genes")

left, right = st.columns(2)


with left:

    top_a = (
        df[df["Direction"] == "Up in A"]
        .nlargest(12, lfc_col)
        .copy()
    )

    fig_a = go.Figure()

    fig_a.add_trace(
        go.Scatter(
            x=top_a[lfc_col],
            y=top_a[gene_col],
            mode="markers",
            marker=dict(
                size=16,
                color="#6D4C91",
                line=dict(
                    width=2,
                    color="white"
                )
            ),
            hovertemplate=(
                "<b>%{y}</b><br>"
                "logFC: %{x:.2f}"
                "<extra></extra>"
            )
        )
    )

    fig_a.update_layout(
        title="Genes enriched in State A",
        height=480,
        xaxis_title="Log fold change",
        yaxis_title="",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
    )

    st.plotly_chart(
        fig_a,
        width="stretch"
    )


with right:

    top_b = (
        df[df["Direction"] == "Up in B"]
        .nsmallest(12, lfc_col)
        .copy()
    )

    top_b["Magnitude"] = (
        top_b[lfc_col].abs()
    )

    fig_b = go.Figure()

    fig_b.add_trace(
        go.Scatter(
            x=top_b["Magnitude"],
            y=top_b[gene_col],
            mode="markers",
            marker=dict(
                size=16,
                color="#D77FA1",
                line=dict(
                    width=2,
                    color="white"
                )
            ),
            hovertemplate=(
                "<b>%{y}</b><br>"
                "|logFC|: %{x:.2f}"
                "<extra></extra>"
            )
        )
    )

    fig_b.update_layout(
        title="Genes enriched in State B",
        height=480,
        xaxis_title="Absolute log fold change",
        yaxis_title="",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
    )

    st.plotly_chart(
        fig_b,
        width="stretch"
    )


# --------------------------------------------------
# SIGNIFICANCE COMPOSITION
# --------------------------------------------------

st.markdown("## DE Signal Composition")

direction_counts = (
    df["Direction"]
    .value_counts()
    .reset_index()
)

direction_counts.columns = [
    "Direction",
    "Genes"
]

fig_donut = px.pie(
    direction_counts,
    names="Direction",
    values="Genes",
    hole=0.62,
    color="Direction",
    color_discrete_map=COLORS,
)

fig_donut.update_traces(
    textinfo="percent+label"
)

fig_donut.update_layout(
    height=420,
    paper_bgcolor="rgba(0,0,0,0)",
    legend_title="",
)

st.plotly_chart(
    fig_donut,
    width="stretch"
)


# --------------------------------------------------
# TABLE
# --------------------------------------------------

with st.expander(
    "Explore complete differential-expression table"
):

    display = df.sort_values(
        padj_col
    )

    st.dataframe(
        display,
        width="stretch",
        hide_index=True,
    )


# --------------------------------------------------
# INTERPRETATION
# --------------------------------------------------

st.markdown("## Interpretation Framework")

c1, c2, c3 = st.columns(3)

with c1:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">AR-HIGH PROGRAM</div>
<div class="sample-state">Prostate lineage</div>
<div class="metric-subtitle">
AR-high comparisons recover characteristic prostate epithelial
genes including KLK and androgen-receptor-associated programs.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c2:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">NEUROENDOCRINE PROGRAM</div>
<div class="sample-state">Neuronal-like transcription</div>
<div class="metric-subtitle">
Neuroendocrine-like states show transcriptional separation
consistent with neuronal and neuroendocrine-associated programs.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c3:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">STEROIDOGENIC PROGRAM</div>
<div class="sample-state">Metabolic specialization</div>
<div class="metric-subtitle">
Steroidogenic-like cells exhibit a distinct transcriptional
program associated with steroid and cholesterol biology.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


st.warning(
    """
    **Important statistical limitation:** these are exploratory cell-state
    differential-expression comparisons. Several states are strongly
    sample-associated, and individual cells are not independent biological
    replicates. Results must not be presented as patient-level differential
    expression.
    """
)

st.caption(
    "CellFlowX • Differential Expression • Wilcoxon • Exploratory cell-state analysis"
)
