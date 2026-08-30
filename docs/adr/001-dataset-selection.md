# ADR-001: Selection of GSE292074 for CellFlowX

## Status

**Accepted — Frozen**

## Context

CellFlowX requires a human single-cell RNA-seq dataset suitable for an end-to-end portfolio project while remaining feasible on local hardware.

The dataset must support:

- genuine single-cell transcriptomic analysis
- biologically meaningful interpretation
- cell-type identification
- tumor heterogeneity analysis
- differential transcriptional analysis
- pathway analysis
- reproducible workflow engineering

The complete CellFlowX project has a target storage ceiling of approximately **1.5 GB**.

## Decision

CellFlowX will use:

**GSE292074 / PRJNA1236646**

The selected study contains human metastatic prostate cancer single-cell RNA-seq data.

Three scRNA-seq samples will be included:

| Sample | GEO Accession | SRA Run |
|---|---|---|
| MH_07-042-M2 | GSM8848584 | SRR32708522 |
| MH_13-084-D12 | GSM8848585 | SRR32708521 |
| MH_13-084-D13 | GSM8848586 | SRR32708520 |

All three selected biological samples will be analyzed.

## Why This Dataset Was Selected

The dataset satisfies the major CellFlowX requirements:

1. Homo sapiens
2. True single-cell RNA sequencing
3. Metastatic prostate cancer biology
4. 10x Genomics data
5. Author-provided processed count matrices
6. Small processed-data footprint
7. Suitable for local Scanpy analysis
8. Supports investigation of inter-tumor heterogeneity
9. Provides a clinically relevant cancer biology story
10. Fits the CellFlowX portfolio objective

## Input Decision

CellFlowX will use the processed Cell Ranger filtered feature-barcode matrices.

Approximate total processed download:

**68.3 MB**

The raw SRA sequencing runs will remain provenance references only.

## Alternatives Considered

Other candidate datasets included:

- gastric cancer scRNA-seq
- breast cancer scRNA-seq
- colorectal cancer scRNA-seq
- lung adenocarcinoma scRNA-seq
- head and neck cancer scRNA-seq
- type 2 diabetes vascular scRNA-seq
- CRISPR single-cell datasets

Some alternatives provided larger cohorts or cleaner case-control comparisons.

However, several had substantially larger data footprints, mixed modalities, weaker comparisons, or greater computational requirements.

GSE292074 provides a strong balance between:

**biological relevance + true scRNA-seq + cancer heterogeneity + computational feasibility**

## Consequences

### Positive

- Complete cohort can be analyzed locally.
- No arbitrary SRR selection is required.
- No raw FASTQ processing is required.
- Storage requirements remain manageable.
- Analysis can focus on advanced downstream single-cell biology.
- Metastatic cancer heterogeneity provides a strong portfolio narrative.

### Limitations

- Only three scRNA-seq tumor samples are available.
- No normal prostate control is included.
- No clean tumor-versus-normal comparison can be performed.
- Sample-level statistical power is limited.
- Population-level conclusions must therefore be conservative.

## Scientific Interpretation

CellFlowX will primarily investigate:

- cellular composition
- inter-tumor heterogeneity
- malignant-cell states
- tumor microenvironment composition
- marker genes
- transcriptional programs
- pathway activity

The project will not present the three samples as a large biological cohort.

## Change Control

This dataset decision is frozen.

Changing the dataset after analysis begins requires:

1. a clear technical or biological justification
2. a new ADR documenting the change
3. an updated sample manifest
4. an updated storage assessment

Dataset changes must not be made silently.
