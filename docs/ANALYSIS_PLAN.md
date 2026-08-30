# CellFlowX Analysis Plan

## Objective

Characterize cellular composition, transcriptional states, and
inter-tumor heterogeneity across metastatic prostate cancer samples
using single-cell RNA sequencing.

## Analysis Workflow

Processed 10x matrices
→ Input validation
→ AnnData construction
→ Cell-level QC
→ Cell/gene filtering
→ Doublet detection
→ Normalization
→ Highly variable genes
→ PCA
→ Neighbor graph
→ Leiden clustering
→ UMAP
→ Marker identification
→ Cell-type annotation
→ Cell composition analysis
→ Malignant-cell analysis
→ Differential analysis
→ Pathway enrichment
→ Biological interpretation

## Stage 1 — Input Validation

Validate:

- expected samples
- matrix/features/barcodes files
- matrix dimensions
- unique cell barcodes
- gene identifiers
- sparse matrix integrity
- metadata consistency

## Stage 2 — QC

Calculate:

- total counts per cell
- genes detected per cell
- mitochondrial percentage
- ribosomal percentage
- doublet score

Filtering thresholds will be determined from observed distributions
rather than blindly hard-coded before data inspection.

## Stage 3 — Preprocessing

Perform:

- normalization
- log transformation
- highly variable gene selection
- scaling where appropriate

Raw counts will be preserved for analyses requiring count data.

## Stage 4 — Dimensionality Reduction

Perform:

- PCA
- neighborhood graph construction
- Leiden clustering
- UMAP

Sample distribution will be inspected before deciding whether batch
correction/integration is justified.

## Stage 5 — Cell-Type Annotation

Clusters will be annotated using:

- cluster marker genes
- established lineage markers
- biological knowledge
- supporting literature where required

Automated annotation will not be accepted without validation.

## Stage 6 — Composition Analysis

Compare:

- cell-type abundance
- cluster abundance
- sample-enriched populations
- sample-depleted populations

across GSM8848584, GSM8848585 and GSM8848586.

## Stage 7 — Malignant-Cell Heterogeneity

Where malignant/epithelial cells can be confidently identified:

- subset malignant cells
- recluster where justified
- identify malignant-cell states
- identify state-specific markers
- compare state distributions across samples
- characterize associated transcriptional programs

## Stage 8 — Differential Analysis

Perform:

- cluster marker analysis
- cell-state comparisons
- sample-associated comparisons where statistically defensible

Cells will not automatically be treated as independent biological
replicates.

Sample-aware or pseudobulk approaches will be considered where
appropriate.

## Stage 9 — Pathway Analysis

Use supported gene signatures for:

- Gene Ontology Biological Process
- Reactome
- MSigDB Hallmark pathways

Pathway interpretation will be driven by observed results.

## Stage 10 — Reporting

Generate:

- QC figures
- UMAPs
- marker plots
- composition plots
- differential-expression summaries
- enrichment plots
- biological conclusions
- limitations

## Validation Principle

Every major stage must pass validation before the workflow proceeds
to the next stage.
