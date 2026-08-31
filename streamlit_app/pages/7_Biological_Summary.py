import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from utils.styles import apply_global_style, metric_card
from utils.paths import HETEROGENEITY_DIR, CNV_DIR, ENRICHMENT_DIR
from utils.loaders import load_csv
from utils.database import (
    get_composition,
    get_heterogeneity,
    get_cnv_summary,
    get_enrichment_summary,
)


st.set_page_config(
    page_title="Biological Summary | CellFlowX",
    page_icon="🧬",
    layout="wide",
)

apply_global_style()

PALETTE = [
    "#6D4C91", "#D77FA1", "#6FA89A",
    "#E89B6D", "#D4AD5C", "#75A6C9",
    "#9C6ADE", "#78A878"
]


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown(
    """
<div class="hero">
<div class="hero-eyebrow">CELLFLOWX • BIOLOGICAL SYNTHESIS</div>
<div class="hero-title">The Tumor Ecosystem at a Glance</div>
<div class="hero-subtitle">
An integrated view connecting cellular composition,
tumor-associated transcriptional states, CNV-like evidence
and functional programs across metastatic prostate cancer.
</div>
</div>
""",
    unsafe_allow_html=True,
)


# --------------------------------------------------
# LOAD
# --------------------------------------------------

heterogeneity = get_heterogeneity().rename(columns={"sample_id": "sample"})

composition_long = get_composition()

composition = composition_long.pivot(
    index="sample_id",
    columns="cell_state",
    values="percentage"
).reset_index().rename(columns={"sample_id": "sample"})

cnv = get_cnv_summary().rename(
    columns={
        "sample_id": "geo_accession",
        "cell_state": "cell_type"
    }
)

enrichment = get_enrichment_summary()

sample_col = (
    "sample"
    if "sample" in composition.columns
    else "geo_accession"
)


# --------------------------------------------------
# TOP KPIs
# --------------------------------------------------

total_cells = int(
    heterogeneity["total_cells"].sum()
)

n_states = len([
    c for c in composition.columns
    if c != sample_col
])

c1, c2, c3, c4 = st.columns(4)

with c1:
    metric_card(
        "TUMORS",
        "3",
        "Metastatic samples"
    )

with c2:
    metric_card(
        "CELLS",
        f"{total_cells:,}",
        "QC-passed atlas"
    )

with c3:
    metric_card(
        "CELL STATES",
        n_states,
        "Expression-defined populations"
    )

with c4:
    metric_card(
        "CORE FINDING",
        "Heterogeneity",
        "Distinct tumor ecosystems"
    )


# --------------------------------------------------
# CENTRAL STORY
# --------------------------------------------------

st.markdown("## Tumor Ecosystem Explorer")

selected_tumor = st.selectbox(
    "Select tumor",
    heterogeneity["sample"].tolist()
)

tumor_summary = heterogeneity[
    heterogeneity["sample"] == selected_tumor
].iloc[0]

tumor_composition = composition[
    composition[sample_col] == selected_tumor
].iloc[0]

state_cols = [
    c for c in composition.columns
    if c != sample_col
]

state_data = pd.DataFrame({
    "Cell state": state_cols,
    "Percentage": [
        float(tumor_composition[c])
        for c in state_cols
    ]
})

state_data = state_data[
    state_data["Percentage"] > 0
].sort_values(
    "Percentage",
    ascending=False
)


left, right = st.columns([1.15, 1])


# --------------------------------------------------
# SUNBURST-LIKE DONUT
# --------------------------------------------------

with left:

    fig = px.pie(
        state_data,
        names="Cell state",
        values="Percentage",
        hole=0.52,
        color_discrete_sequence=PALETTE,
        title=f"{selected_tumor} cellular ecosystem",
    )

    fig.update_traces(
        textinfo="percent",
        hovertemplate=(
            "<b>%{label}</b><br>"
            "%{value:.2f}%"
            "<extra></extra>"
        )
    )

    fig.update_layout(
        height=520,
        paper_bgcolor="rgba(0,0,0,0)",
        legend_title="Cell state"
    )

    st.plotly_chart(
        fig,
        width="stretch"
    )


with right:

    st.markdown(
        f"""
<div class="glass-card">
<div class="card-title">Tumor identity</div>

<div class="sample-id">{selected_tumor}</div>

<div class="sample-state">
{tumor_summary['dominant_state']}
</div>

<div class="sample-pct">
{tumor_summary['dominant_state_pct']:.1f}%
</div>

<div class="metric-subtitle">
dominant cellular state
</div>

<br>

<div class="card-text">
<b>{int(tumor_summary['total_cells']):,}</b> cells were retained
from this tumor, representing
<b>{int(tumor_summary['detected_states'])}</b> detected cellular states.

The composition profile highlights substantial differences
between metastatic tumor ecosystems.
</div>
</div>
""",
        unsafe_allow_html=True,
    )


# --------------------------------------------------
# BIOLOGICAL AXES
# --------------------------------------------------

st.markdown("## Four Biological Axes")

a1, a2, a3, a4 = st.columns(4)

with a1:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">01 • CELLULAR ECOLOGY</div>
<div class="sample-state">Composition</div>
<div class="metric-subtitle">
Tumors differ substantially in epithelial, neuroendocrine,
stromal and immune composition.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with a2:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">02 • CELL STATE</div>
<div class="sample-state">Transcription</div>
<div class="metric-subtitle">
AR-high, luminal, neuroendocrine-like and steroidogenic-like
programs define distinct tumor-associated states.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with a3:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">03 • CNV EVIDENCE</div>
<div class="sample-state">Genome proxy</div>
<div class="metric-subtitle">
Chromosome-level expression deviations provide supportive
evidence for selected epithelial populations.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with a4:
    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">04 • FUNCTION</div>
<div class="sample-state">Pathways</div>
<div class="metric-subtitle">
Functional enrichment links transcriptional states to
distinct biological programs.
</div>
</div>
""",
        unsafe_allow_html=True
    )


# --------------------------------------------------
# CNV EVIDENCE MATRIX
# --------------------------------------------------

st.markdown("## Tumor-State Evidence Matrix")

required = {
    "geo_accession",
    "cell_type",
    "cnv_high_pct"
}

if required.issubset(cnv.columns):

    cnv_heat = cnv.pivot_table(
        index="cell_type",
        columns="geo_accession",
        values="cnv_high_pct",
        fill_value=0
    )

    fig_cnv = px.imshow(
        cnv_heat,
        aspect="auto",
        text_auto=".1f",
        color_continuous_scale=[
            "#FBF8FD",
            "#E2D2EC",
            "#B895CE",
            "#8056A0",
            "#4D2D63"
        ],
        title="CNV-like expression evidence across tumor-associated states"
    )

    fig_cnv.update_layout(
        height=520,
        xaxis_title="Tumor",
        yaxis_title="Cell state",
        coloraxis_colorbar_title="CNV-high %",
        paper_bgcolor="rgba(0,0,0,0)"
    )

    st.plotly_chart(
        fig_cnv,
        width="stretch"
    )


# --------------------------------------------------
# FUNCTIONAL PROGRAM NETWORK-LIKE VIEW
# --------------------------------------------------

st.markdown("## Functional Program Map")

functional = enrichment.copy()

functional["Program"] = (
    functional["gene_set"]
    .str.replace("_up_A", " ↑A", regex=False)
    .str.replace("_up_B", " ↑B", regex=False)
    .str.replace("_", " ", regex=False)
)

functional["Significance"] = (
    -functional["top_adjusted_p"]
    .clip(lower=1e-300)
    .apply(lambda x: __import__("math").log10(x))
)

fig_func = px.scatter(
    functional,
    x="input_genes",
    y="Significance",
    size="significant_pathways",
    color="Program",
    hover_name="top_pathway",
    hover_data={
        "input_genes": True,
        "significant_pathways": True,
        "top_adjusted_p": ":.2e",
        "Significance": False
    },
    color_discrete_sequence=PALETTE,
    size_max=65,
    title="Transcriptional programs → functional biology",
)

fig_func.update_layout(
    height=540,
    xaxis_title="Differential genes",
    yaxis_title="Top pathway significance (−log10 FDR)",
    legend_title="Gene program",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
)

st.plotly_chart(
    fig_func,
    width="stretch"
)


# --------------------------------------------------
# MAIN FINDINGS
# --------------------------------------------------

st.markdown("## Key Biological Findings")

c1, c2, c3 = st.columns(3)

with c1:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">TUMOR 1 • GSM8848584</div>
<div class="sample-state">AR-high + steroidogenic</div>
<div class="metric-subtitle">
Steroidogenic-like cells represent 32.7% and AR-high prostate
epithelial cells 32.2%, creating a strongly epithelial/metabolic ecosystem.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with c2:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">TUMOR 2 • GSM8848585</div>
<div class="sample-state">Neuroendocrine-rich</div>
<div class="metric-subtitle">
Neuroendocrine-like cells dominate at 43.5%, accompanied by
luminal epithelial and cycling populations.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with c3:

    st.markdown(
        """
<div class="sample-card">
<div class="sample-id">TUMOR 3 • GSM8848586</div>
<div class="sample-state">Stromal / myeloid-rich</div>
<div class="metric-subtitle">
Activated fibroblasts dominate at 56.1%, with macrophages
forming another major component of the tumor ecosystem.
</div>
</div>
""",
        unsafe_allow_html=True
    )


# --------------------------------------------------
# MOLECULAR PROGRAMS
# --------------------------------------------------

st.markdown("## Molecular Interpretation")

m1, m2, m3 = st.columns(3)

with m1:

    st.markdown(
        """
<div class="glass-card">
<div class="card-title">Neuroendocrine-like</div>
<div class="card-text">
Marker and differential-expression analyses identify neuronal
and neuroendocrine-associated transcription, while enrichment highlights
<b>nervous system development</b> and
<b>chemical synaptic transmission</b>.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with m2:

    st.markdown(
        """
<div class="glass-card">
<div class="card-title">Steroidogenic-like</div>
<div class="card-text">
Genes including steroidogenic machinery define this state,
with functional enrichment strongly highlighting
<b>cholesterol metabolic processes</b>.
</div>
</div>
""",
        unsafe_allow_html=True
    )

with m3:

    st.markdown(
        """
<div class="glass-card">
<div class="card-title">Prostate epithelial</div>
<div class="card-text">
AR-high and luminal states retain prostate-lineage transcriptional
programs. CNV-like expression deviations provide additional supportive
evidence particularly within these epithelial populations.
</div>
</div>
""",
        unsafe_allow_html=True
    )


# --------------------------------------------------
# FINAL MODEL
# --------------------------------------------------

st.markdown("## CellFlowX Biological Model")

st.markdown(
    """
<div class="glass-card">
<div class="card-title">Integrated interpretation</div>
<div class="card-text">

<b>Metastatic prostate cancer in this dataset does not present as one
uniform cellular state.</b>

Instead, the three tumors occupy markedly different cellular ecosystems:

<br><br>

<b>GSM8848584</b> → AR-high / steroidogenic-rich<br>
<b>GSM8848585</b> → neuroendocrine / luminal-rich<br>
<b>GSM8848586</b> → activated fibroblast / macrophage-rich

<br><br>

These compositional differences are accompanied by distinct
transcriptional programs and functional enrichment patterns.
Expression-based CNV analysis provides additional supportive evidence
for genomic dysregulation in selected prostate epithelial populations.

</div>
</div>
""",
    unsafe_allow_html=True
)


# --------------------------------------------------
# LIMITATIONS
# --------------------------------------------------

with st.expander(
    "Scientific interpretation & limitations",
    expanded=False
):

    st.markdown(
        """
        **Important boundaries of interpretation**

        - The study contains only **three metastatic tumor samples**.
        - There is **no normal/control cohort**.
        - Several cell states are strongly associated with individual samples.
        - Cells are not independent biological replicates.
        - Differential expression is therefore exploratory cell-state analysis.
        - Pathway enrichment inherits this sample-confounding limitation.
        - CNV scores are chromosome-level **RNA-expression proxies**, not DNA CNV calls.
        - Expression-defined states do not by themselves prove malignant identity.
        - UMAP represents transcriptional similarity, not spatial or lineage relationships.

        CellFlowX therefore demonstrates **within-dataset biological
        heterogeneity and reproducible single-cell analysis**, rather than
        population-level clinical inference.
        """
    )


st.caption(
    "CellFlowX • Biological Synthesis • Metastatic Prostate Cancer • Single-Cell Transcriptomics"
)
