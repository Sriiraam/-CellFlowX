import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results" / "provenance"
OUTDIR.mkdir(parents=True, exist_ok=True)


def command(cmd):
    try:
        return subprocess.check_output(
            cmd,
            shell=True,
            text=True,
            stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unavailable"


payload = {
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "project": "CellFlowX",
    "python_version": sys.version.split()[0],
    "platform": platform.platform(),
    "git_commit": command("git rev-parse HEAD"),
    "git_branch": command("git rev-parse --abbrev-ref HEAD"),
    "nextflow_version": command("nextflow -version | grep -m1 version"),
    "java_version": command("java -version 2>&1 | head -n 1"),
    "docker_version": command("docker --version"),
    "working_directory": str(ROOT),
    "dataset": {
        "geo_series": "GSE292074",
        "bioproject": "PRJNA1236646",
        "samples": [
            "GSM8848584",
            "GSM8848585",
            "GSM8848586"
        ]
    },
    "analysis": {
        "raw_cells": 8660,
        "qc_retained_cells": 8233,
        "hvg_count": 2000,
        "leiden_clusters": 15
    }
}

with open(OUTDIR / "provenance.json", "w") as handle:
    json.dump(payload, handle, indent=2)

requirements = ROOT / "requirements.txt"

if requirements.exists():
    (OUTDIR / "requirements_snapshot.txt").write_text(
        requirements.read_text()
    )

config = ROOT / "nextflow.config"

if config.exists():
    (OUTDIR / "nextflow_config_snapshot.txt").write_text(
        config.read_text()
    )

print(f"Provenance written to: {OUTDIR}")
