from pathlib import Path
import pandas as pd
import re

ROOT = Path(__file__).resolve().parents[1]

TRACE = ROOT / "results" / "benchmarks" / "nextflow_trace.txt"
OUTDIR = ROOT / "results" / "benchmarks"
OUTDIR.mkdir(parents=True, exist_ok=True)


def parse_memory(value):
    if pd.isna(value):
        return None

    text = str(value).strip()

    units = {
        "B": 1 / (1024 ** 3),
        "KB": 1 / (1024 ** 2),
        "MB": 1 / 1024,
        "GB": 1,
        "TB": 1024,
    }

    m = re.match(r"([\d.]+)\s*([KMGT]?B)", text, re.I)

    if not m:
        return None

    number = float(m.group(1))
    unit = m.group(2).upper()

    return number * units[unit]


def parse_duration(value):
    if pd.isna(value):
        return None

    text = str(value).strip()

    total = 0.0

    patterns = {
        "d": 86400,
        "h": 3600,
        "m": 60,
        "s": 1,
        "ms": 0.001,
    }

    for number, unit in re.findall(
        r"([\d.]+)\s*(ms|d|h|m|s)",
        text
    ):
        total += float(number) * patterns[unit]

    return total


if not TRACE.exists():
    raise SystemExit(
        f"Trace file not found: {TRACE}"
    )

df = pd.read_csv(
    TRACE,
    sep="\t"
)

print("Trace columns:")
print(", ".join(df.columns))


# --------------------------------------------------
# NORMALISE PROCESS NAME
# --------------------------------------------------

process_col = (
    "process"
    if "process" in df.columns
    else "name"
    if "name" in df.columns
    else None
)

if process_col is None:
    raise SystemExit(
        "Nextflow trace contains neither 'process' nor 'name'."
    )

df["process_clean"] = (
    df[process_col]
    .astype(str)
    .str.replace(r"\s*\(.*\)$", "", regex=True)
)


# --------------------------------------------------
# DURATION
# --------------------------------------------------

duration_col = None

for candidate in [
    "realtime",
    "duration"
]:
    if candidate in df.columns:
        duration_col = candidate
        break

if duration_col:
    df["runtime_seconds"] = df[
        duration_col
    ].apply(parse_duration)
else:
    df["runtime_seconds"] = None


# --------------------------------------------------
# MEMORY
# --------------------------------------------------

memory_col = None

for candidate in [
    "peak_rss",
    "rss"
]:
    if candidate in df.columns:
        memory_col = candidate
        break

if memory_col:
    df["peak_memory_gb"] = df[
        memory_col
    ].apply(parse_memory)
else:
    df["peak_memory_gb"] = None


# --------------------------------------------------
# CPU
# --------------------------------------------------

if "%cpu" in df.columns:
    df["cpu_percent"] = (
        df["%cpu"]
        .astype(str)
        .str.replace("%", "", regex=False)
    )

    df["cpu_percent"] = pd.to_numeric(
        df["cpu_percent"],
        errors="coerce"
    )
else:
    df["cpu_percent"] = None


# --------------------------------------------------
# PROCESS SUMMARY
# --------------------------------------------------

summary = (
    df.groupby(
        "process_clean",
        dropna=False
    )
    .agg(
        tasks=("process_clean", "size"),
        runtime_seconds=("runtime_seconds", "sum"),
        mean_runtime_seconds=("runtime_seconds", "mean"),
        peak_memory_gb=("peak_memory_gb", "max"),
        mean_cpu_percent=("cpu_percent", "mean"),
    )
    .reset_index()
)

summary = summary.sort_values(
    "runtime_seconds",
    ascending=False
)

summary.to_csv(
    OUTDIR / "process_benchmark_summary.csv",
    index=False
)


# --------------------------------------------------
# PROJECT SUMMARY
# --------------------------------------------------

total_runtime = df[
    "runtime_seconds"
].sum()

peak_memory = df[
    "peak_memory_gb"
].max()

mean_cpu = df[
    "cpu_percent"
].mean()

project = pd.DataFrame([
    {
        "pipeline": "CellFlowX",
        "tasks_recorded": len(df),
        "processes_recorded": df["process_clean"].nunique(),
        "total_task_runtime_seconds": round(total_runtime, 2),
        "total_task_runtime_minutes": round(total_runtime / 60, 2),
        "peak_task_memory_gb": (
            round(peak_memory, 3)
            if pd.notna(peak_memory)
            else None
        ),
        "mean_task_cpu_percent": (
            round(mean_cpu, 2)
            if pd.notna(mean_cpu)
            else None
        ),
    }
])

project.to_csv(
    OUTDIR / "project_benchmark_summary.csv",
    index=False
)


# --------------------------------------------------
# MARKDOWN
# --------------------------------------------------

md = [
    "# CellFlowX Benchmark Summary",
    "",
    "Benchmark metrics were derived from the Nextflow execution trace.",
    "",
    "## Pipeline Summary",
    "",
    f"- Tasks recorded: {len(df)}",
    f"- Processes recorded: {df['process_clean'].nunique()}",
    f"- Total task runtime: {total_runtime / 60:.2f} minutes",
]

if pd.notna(peak_memory):
    md.append(
        f"- Peak task memory: {peak_memory:.2f} GB"
    )

if pd.notna(mean_cpu):
    md.append(
        f"- Mean task CPU utilisation: {mean_cpu:.1f}%"
    )

md += [
    "",
    "## Process-Level Performance",
    "",
    summary.to_markdown(index=False),
    "",
    "### Interpretation",
    "",
    "These metrics represent process-level resource usage reported by Nextflow "
    "on the local CellFlowX execution environment. They are intended for "
    "workflow engineering and reproducibility assessment rather than "
    "cross-platform performance claims.",
]

(
    OUTDIR /
    "BENCHMARKS.md"
).write_text(
    "\n".join(md)
)

print()
print("✓ process_benchmark_summary.csv")
print("✓ project_benchmark_summary.csv")
print("✓ BENCHMARKS.md")
print()
print(project.to_string(index=False))
