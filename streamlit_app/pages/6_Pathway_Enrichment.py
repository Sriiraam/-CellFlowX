import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

from utils.styles import apply_global_style, metric_card
from utils.paths import ENRICHMENT_DIR
from utils.loaders import load_csv


st.set_page_config(
    page_title="Pathway Enrichment | CellFlowX",
    page_icon="🧬",
    layout="wide",
)

apply_global_style()

PALETTE = [
    "#6D4C91", "#9C6ADE", "#D77FA1", "#E89B6D",
    "#D4AD5C", "#6FA89A", "#75A6C9", "#A98CC1",
]


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">CELLFLOWX • FUNCTIONAL ENRICHMENT</div>
<div class="hero-title">From Genes to Biological Programs</div>
<div class="hero-subtitle">
Explore GO Biological Process enrichment underlying AR-high,
luminal, neuroendocrine-like and steroidogenic-like transcriptional states.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# LOAD SUMMARY
# --------------------------------------------------

summary_path = (
    ENRICHMENT_DIR /
    "functional_enrichment_summary.csv"
)

summary = load_csv(summary_path)

summary["Display"] = (
    summary["gene_set"]
    .str.replace("_up_A", " • State A", regex=False)
    .str.replace("_up_B", " • State B", regex=False)
    .str.replace("_", " ", regex=False)
)


# --------------------------------------------------
# KPIs
# --------------------------------------------------

total_programs = len(summary)

total_pathways = int(
    summary["significant_pathways"].sum()
)

strongest = summary.loc[
    summary["top_adjusted_p"].idxmin()
]

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "GENE SETS",
        total_programs,
        "Directional DE programs"
    )

with c2:
    metric_card(
        "GO TERMS",
        f"{total_pathways:,}",
        "Significant enrichments"
    )

with c3:
    metric_card(
        "DATABASE",
        "GO BP",
        "Biological Process 2023"
    )

with c4:
    metric_card(
        "STRONGEST FDR",
        f"{strongest['top_adjusted_p']:.1e}",
        "Top enrichment signal"
    )


# --------------------------------------------------
# PROGRAM SELECTOR
# --------------------------------------------------

st.markdown("## Functional Program Explorer")

selected_gene_set = st.selectbox(
    "Choose directional gene program",
    summary["gene_set"].tolist(),
    format_func=lambda x: (
        x.replace("_up_A", " • Up in State A")
         .replace("_up_B", " • Up in State B")
         .replace("_", " ")
    )
)


# --------------------------------------------------
# LOCATE ENRICHMENT FILE
# --------------------------------------------------

csv_files = [
    f for f in ENRICHMENT_DIR.glob("*.csv")
    if f.name != "functional_enrichment_summary.csv"
]

matching = [
    f for f in csv_files
    if selected_gene_set in f.stem
]

if not matching:
    st.warning(
        "Could not locate the enrichment table for this program."
    )
    st.stop()

enrich = load_csv(matching[0])


def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


term_col = find_col(
    enrich,
    ["Term", "term", "Pathway", "pathway"]
)

padj_col = find_col(
    enrich,
    [
        "Adjusted P-value",
        "Adjusted P-value ",
        "adjusted_pvalue",
        "padj",
        "FDR",
    ]
)

overlap_col = find_col(
    enrich,
    ["Overlap", "overlap"]
)

combined_col = find_col(
    enrich,
    ["Combined Score", "combined_score"]
)


if term_col is None or padj_col is None:
    st.error(
        "Required enrichment columns could not be detected."
    )
    st.stop()


enrich[padj_col] = pd.to_numeric(
    enrich[padj_col],
    errors="coerce"
)

enrich = enrich.dropna(
    subset=[padj_col]
)

enrich["Significance"] = -np.log10(
    enrich[padj_col].clip(lower=1e-300)
)


# Parse Enrichr overlap such as 25/300
if overlap_col:

    enrich["Genes"] = (
        enrich[overlap_col]
        .astype(str)
        .str.split("/")
        .str[0]
    )

    enrich["Genes"] = pd.to_numeric(
        enrich["Genes"],
        errors="coerce"
    )

else:
    enrich["Genes"] = 1


# --------------------------------------------------
# SELECTED PROGRAM KPIs
# --------------------------------------------------

selected_summary = summary[
    summary["gene_set"] == selected_gene_set
].iloc[0]

s1, s2, s3 = st.columns(3)

s1.metric(
    "Input DE genes",
    f"{int(selected_summary['input_genes']):,}"
)

s2.metric(
    "Significant GO terms",
    int(selected_summary["significant_pathways"])
)

s3.metric(
    "Top program",
    selected_summary["top_pathway"]
)


# --------------------------------------------------
# BUBBLE LANDSCAPE
# --------------------------------------------------

st.markdown("## Pathway Landscape")

top = (
    enrich
    .sort_values(padj_col)
    .head(20)
    .copy()
)

if combined_col:
    top["Bubble"] = pd.to_numeric(
        top[combined_col],
        errors="coerce"
    ).abs()

    top["Bubble"] = top["Bubble"].fillna(
        top["Genes"]
    )
else:
    top["Bubble"] = top["Genes"]


fig_bubble = px.scatter(
    top,
    x="Significance",
    y=term_col,
    size="Bubble",
    color="Significance",
    hover_data={
        padj_col: ":.2e",
        "Genes": True,
        "Bubble": False,
    },
    color_continuous_scale=[
        "#E8DDF1",
        "#C39DDB",
        "#9568B6",
        "#5B376F",
    ],
    size_max=42,
    title="Top enriched biological processes",
)

fig_bubble.update_layout(
    height=680,
    xaxis_title="Enrichment significance (−log10 adjusted P)",
    yaxis_title="",
    coloraxis_colorbar_title="Significance",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(
    fig_bubble,
    width="stretch"
)


# --------------------------------------------------
# GLOBAL PROGRAM LANDSCAPE
# --------------------------------------------------

st.markdown("## Cross-State Functional Landscape")

global_df = summary.copy()

global_df["Significance"] = -np.log10(
    global_df["top_adjusted_p"].clip(lower=1e-300)
)

fig_global = px.scatter(
    global_df,
    x="input_genes",
    y="Significance",
    size="significant_pathways",
    color="gene_set",
    hover_name="top_pathway",
    hover_data={
        "input_genes": True,
        "significant_pathways": True,
        "top_adjusted_p": ":.2e",
        "Significance": False,
    },
    color_discrete_sequence=PALETTE,
    size_max=55,
    title="Functional complexity across tumor-state programs",
)

fig_global.update_layout(
    height=530,
    xaxis_title="Differentially expressed genes entering enrichment",
    yaxis_title="Top pathway significance (−log10 FDR)",
    legend_title="Gene program",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(
    fig_global,
    width="stretch"
)


# --------------------------------------------------
# TOP PATHWAY CARDS
# --------------------------------------------------

st.markdown("## Major Biological Programs")

c1, c2, c3 = st.columns(3)

with c1:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">NEUROENDOCRINE-LIKE</div>
<div class="sample-state">Neural programs</div>
<div class="metric-subtitle">
Nervous system development and chemical synaptic transmission
are strongly enriched in neuroendocrine-associated transcription.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c2:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">STEROIDOGENIC-LIKE</div>
<div class="sample-state">Cholesterol metabolism</div>
<div class="metric-subtitle">
Steroidogenic-like transcription is strongly associated with
cholesterol metabolic processes.
</div>
</div>
""",
        unsafe_allow_html=True,
    )

with c3:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">LUMINAL PROGRAM</div>
<div class="sample-state">Cell migration biology</div>
<div class="metric-subtitle">
Luminal-enriched genes show functional association with
regulation of cellular migration.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


# --------------------------------------------------
# TOP TERM EXPLORER
# --------------------------------------------------

st.markdown("## Pathway Inspector")

pathway_options = top[term_col].astype(str).tolist()

selected_pathway = st.selectbox(
    "Select GO biological process",
    pathway_options
)

pathway_row = top[
    top[term_col].astype(str) == selected_pathway
].iloc[0]

p1, p2, p3 = st.columns(3)

p1.metric(
    "Adjusted P-value",
    f"{pathway_row[padj_col]:.2e}"
)

p2.metric(
    "−log10 FDR",
    f"{pathway_row['Significance']:.2f}"
)

p3.metric(
    "Overlapping genes",
    int(pathway_row["Genes"])
    if pd.notna(pathway_row["Genes"])
    else "-"
)


# --------------------------------------------------
# TABLE
# --------------------------------------------------

with st.expander(
    "Explore complete enrichment results"
):

    st.dataframe(
        enrich.sort_values(padj_col),
        width="stretch",
        hide_index=True,
    )


# --------------------------------------------------
# LIMITATION
# --------------------------------------------------

st.warning(
    """
    **Interpretation limitation:** pathway enrichment inherits the limitations
    of the Phase 10 exploratory cell-state differential-expression analysis.
    Several cellular states are strongly sample-associated, so enrichment
    describes transcriptional programs observed in this dataset rather than
    patient-independent pathway effects.
    """
)

st.caption(
    "CellFlowX • Functional Enrichment • GO Biological Process • GSEApy / Enrichr"
)
