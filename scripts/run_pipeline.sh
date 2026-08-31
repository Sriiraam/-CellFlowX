#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs results/provenance

STAMP=$(date +"%Y%m%d_%H%M%S")
LOG="logs/cellflowx_${STAMP}.log"

echo "CellFlowX run started: $(date)" | tee "$LOG"

nextflow run main.nf \
    -resume \
    "$@" 2>&1 | tee -a "$LOG"

python scripts/record_provenance.py 2>&1 | tee -a "$LOG"

echo "CellFlowX run completed: $(date)" | tee -a "$LOG"
echo "Log: $LOG"
