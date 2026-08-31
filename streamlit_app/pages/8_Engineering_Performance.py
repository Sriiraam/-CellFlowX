import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "streamlit_app"))

from utils.styles import apply_global_style

apply_global_style()

st.set_page_config(
    page_title="Engineering & Performance | CellFlowX",
    page_icon="⚙️",
    layout="wide",
)

BENCH = ROOT / "streamlit_app" / "data" / "results" / "benchmarks"

project = pd.read_csv(BENCH / "project_benchmark_summary.csv")
processes = pd.read_csv(BENCH / "process_benchmark_summary.csv")

row = project.iloc[0]

st.markdown("""
<div style="
    padding: 2rem 2.2rem;
    border-radius: 24px;
    background: linear-gradient(135deg,#241238,#51306f,#76539a);
    color: white;
    margin-bottom: 1.6rem;
    box-shadow: 0 12px 32px rgba(62,35,85,0.18);
">
    <div style="font-size:0.85rem;letter-spacing:0.12em;
                text-transform:uppercase;opacity:0.8;">
        Workflow Engineering
    </div>
    <div style="font-size:2.15rem;font-weight:750;margin-top:0.3rem;">
        Engineering & Performance
    </div>
    <div style="font-size:1.02rem;opacity:0.88;margin-top:0.45rem;">
        Reproducibility, resource utilisation and process-level
        performance of the CellFlowX workflow.
    </div>
</div>
""", unsafe_allow_html=True)

c1, c2, c3, c4 = st.columns(4)

c1.metric("Processes", int(row["processes_recorded"]))
c2.metric("Tasks", int(row["tasks_recorded"]))
c3.metric(
    "Summed Task Runtime",
    f'{row["total_task_runtime_minutes"]:.2f} min'
)
c4.metric(
    "Peak Task Memory",
    f'{row["peak_task_memory_gb"]:.2f} GB'
)

st.caption(
    "Summed task runtime represents cumulative process execution time "
    "from the Nextflow trace and should not be interpreted as pipeline "
    "wall-clock runtime."
)

st.markdown("### Runtime profile")

runtime = processes.sort_values(
    "runtime_seconds",
    ascending=True
)

fig_runtime = px.bar(
    runtime,
    x="runtime_seconds",
    y="process_clean",
    orientation="h",
    color="runtime_seconds",
    color_continuous_scale=[
        "#D8C4EA",
        "#A779C8",
        "#68407F",
    ],
    labels={
        "runtime_seconds": "Runtime (seconds)",
        "process_clean": "Process",
    },
)

fig_runtime.update_layout(
    coloraxis_showscale=False,
    height=560,
    margin=dict(l=20, r=20, t=20, b=20),
)

st.plotly_chart(fig_runtime, use_container_width=True)

st.markdown("### Resource utilisation")

left, right = st.columns(2)

with left:
    memory = processes.sort_values(
        "peak_memory_gb",
        ascending=False
    )

    fig_mem = px.scatter(
        memory,
        x="runtime_seconds",
        y="peak_memory_gb",
        size="peak_memory_gb",
        color="process_clean",
        hover_name="process_clean",
        labels={
            "runtime_seconds": "Runtime (seconds)",
            "peak_memory_gb": "Peak memory (GB)",
            "process_clean": "Process",
        },
    )

    fig_mem.update_layout(
        showlegend=False,
        height=430,
        margin=dict(l=20, r=20, t=20, b=20),
    )

    st.plotly_chart(fig_mem, use_container_width=True)

with right:
    cpu = processes.sort_values(
        "mean_cpu_percent",
        ascending=False
    ).head(8)

    fig_cpu = px.bar(
        cpu,
        x="process_clean",
        y="mean_cpu_percent",
        color="mean_cpu_percent",
        color_continuous_scale=[
            "#F4D6C9",
            "#D88C72",
            "#8A4F66",
        ],
        labels={
            "process_clean": "Process",
            "mean_cpu_percent": "Mean CPU utilisation (%)",
        },
    )

    fig_cpu.update_layout(
        coloraxis_showscale=False,
        height=430,
        xaxis_tickangle=-45,
        margin=dict(l=20, r=20, t=20, b=20),
    )

    st.plotly_chart(fig_cpu, use_container_width=True)

st.markdown("### Process benchmark table")

display = processes[
    [
        "process_clean",
        "runtime_seconds",
        "peak_memory_gb",
        "mean_cpu_percent",
    ]
].copy()

display.columns = [
    "Process",
    "Runtime (s)",
    "Peak Memory (GB)",
    "Mean CPU (%)",
]

display["Runtime (s)"] = display["Runtime (s)"].round(1)
display["Peak Memory (GB)"] = display["Peak Memory (GB)"].round(2)
display["Mean CPU (%)"] = display["Mean CPU (%)"].round(1)

st.dataframe(
    display,
    use_container_width=True,
    hide_index=True,
)

st.markdown("### Engineering interpretation")

st.info(
    "CellFlowX completed all recorded analytical stages within modest "
    "memory requirements on a local WSL2 environment. Functional "
    "enrichment and embedding/clustering were the largest runtime "
    "contributors, while peak task memory remained approximately "
    "1.3 GB. These measurements describe this specific local execution "
    "and are not intended as cross-platform performance claims."
)
