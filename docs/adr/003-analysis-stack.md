# ADR-003: CellFlowX Analysis and Engineering Stack

## Status

**Accepted — Frozen**

## Context

CellFlowX requires a technology stack capable of performing
single-cell RNA-seq analysis while demonstrating reproducible
bioinformatics workflow engineering.

The workflow must remain feasible on local hardware.

## Decision

### Primary Language

**Python**

Python will be the primary analysis and engineering language.

### Single-Cell Framework

**Scanpy + AnnData**

Scanpy will provide the primary downstream single-cell analysis
framework.

AnnData will provide the core single-cell data structure.

### Supporting Libraries

Expected supporting libraries include:

- pandas
- NumPy
- SciPy
- matplotlib
- plotly
- scrublet
- gseapy and/or decoupler

Additional libraries will only be introduced when required.

## Workflow Orchestration

**Nextflow DSL2**

Nextflow will orchestrate reproducible analysis stages.

The workflow will use modular processes and subworkflows where
appropriate.

## Containerization

**Docker**

Docker will provide a reproducible execution environment.

## Testing

**pytest**

Automated tests will validate important workflow and input-handling
functions.

## Continuous Integration

**GitHub Actions**

CI will execute lightweight validation and testing.

Large biological datasets will not be downloaded during routine CI.

## Dashboard

**Streamlit**

Streamlit will provide an interactive biological results dashboard.

## R / Seurat

R and Seurat are not part of the default CellFlowX stack.

They may only be introduced if they provide a clearly justified
analysis unavailable or materially weaker in the Python workflow.

## SQL

SQLite may be introduced for lightweight dashboard/result querying
only if it provides a genuine engineering advantage.

SQL will not be added merely to increase the number of technologies
used by the project.

## HPC and Cloud

CellFlowX will execute locally first.

The workflow architecture may later support:

- SLURM
- HPC
- containerized cloud execution

Paid cloud compute is not required for project completion.

## Decision Principle

Technologies are included only when they improve:

- reproducibility
- maintainability
- biological analysis
- testing
- portability
- usability

Unnecessary technologies will not be added solely for portfolio
decoration.
