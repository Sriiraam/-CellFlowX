
CellFlowX Benchmark Plan
Objective

Measure the computational requirements of CellFlowX on local
hardware.

Metrics

Record:

total runtime
runtime per major stage
peak memory where measurable
input storage
intermediate storage
final project storage
number of cells processed
number of genes
number of cells retained after QC
Nextflow Metrics

Where applicable, retain:

trace
execution report
timeline
Reproducibility Benchmark

A clean rerun should reproduce the expected workflow outputs using
the frozen configuration and inputs.

Resume Benchmark

Nextflow resume behaviour will be tested to demonstrate that
completed stages are not unnecessarily recomputed.

Storage Benchmark

The complete project should remain approximately:

≤ 1.5 GB

Storage will be checked at major milestones.

Purpose

Benchmarking demonstrates that CellFlowX is not only biologically
functional but also computationally characterized and operationally
reproducible.
