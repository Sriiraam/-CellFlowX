# CellFlowX Project Scope

## Project Title

**CellFlowX: Single-Cell Transcriptomic Profiling of Tumor Heterogeneity in Metastatic Prostate Cancer**

## Project Objective

CellFlowX is an end-to-end, reproducible single-cell RNA-seq analysis workflow designed to characterize cellular and transcriptional heterogeneity in human metastatic prostate cancer.

The project combines biological analysis with production-style bioinformatics engineering.

## Dataset

- GEO Series: GSE292074
- BioProject: PRJNA1236646
- Organism: Homo sapiens
- Assay: Single-cell RNA sequencing
- Tissue: Metastatic prostate cancer tumor tissue
- Platform: 10x Genomics
- Input type: Processed Cell Ranger filtered feature-barcode matrices

## Frozen Samples

- GSM8848584 — MH_07-042-M2
- GSM8848585 — MH_13-084-D12
- GSM8848586 — MH_13-084-D13

## Primary Biological Question

How do cellular composition and transcriptional states vary across metastatic prostate cancer samples, and which cell populations and molecular programs contribute to inter-tumor heterogeneity?

## Planned Analysis

1. Input validation
2. Single-cell QC
3. Cell and gene filtering
4. Mitochondrial QC
5. Doublet detection
6. Normalization
7. Highly variable gene selection
8. PCA
9. Neighbor graph construction
10. Leiden clustering
11. UMAP
12. Marker-gene identification
13. Cell-type annotation
14. Cell-population composition analysis
15. Malignant-cell heterogeneity analysis
16. Differential expression
17. Pathway and enrichment analysis
18. Biological interpretation
19. Benchmarking
20. Interactive Streamlit dashboard

## Engineering Scope

CellFlowX will include:

- Python
- Scanpy
- AnnData
- Nextflow DSL2
- Docker
- pytest
- GitHub Actions
- Streamlit
- Git/GitHub
- reproducible configuration
- runtime and storage benchmarking

## Data Strategy

Raw FASTQ/SRA data will not be downloaded.

CellFlowX will use author-provided processed 10x Genomics count matrices.

This decision keeps the complete project within the defined local storage and hardware limits while preserving the biological single-cell analysis workflow.

## Reference Strategy

No reference genome or alignment index is required for the main CellFlowX workflow because read alignment and gene quantification were already performed by the original study using Cell Ranger.

Therefore CellFlowX will not download:

- FASTA
- GTF
- STAR index
- Cell Ranger reference bundle

## Storage Constraint

The complete CellFlowX project should remain approximately:

**≤ 1.5 GB**

Expected processed dataset download:

**~68.3 MB**

## Scope Limitation

The three selected samples are metastatic prostate cancer samples.

CellFlowX will therefore focus on:

- inter-tumor heterogeneity
- cellular composition
- malignant-cell states
- tumor microenvironment differences

It will not claim a tumor-vs-normal comparison.

## Execution Principle

Every major analysis stage must be validated before proceeding to the next stage.

The project will be developed locally first and designed to remain reproducible and portable for future HPC or cloud execution.