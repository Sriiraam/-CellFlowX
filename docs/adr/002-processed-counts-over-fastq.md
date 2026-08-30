# ADR-002: Use Processed 10x Count Matrices Instead of Raw FASTQ

## Status

**Accepted — Frozen**

## Context

CellFlowX is designed as a downstream single-cell RNA-seq analysis and bioinformatics engineering project.

The selected dataset, GSE292074, provides both sequencing provenance through SRA and author-processed 10x Genomics feature-barcode matrices.

Processing the raw sequencing data would require substantially more storage, compute, and additional reference resources.

CellFlowX has a complete project storage target of:

**≤ 1.5 GB**

## Raw Sequencing Data

The corresponding SRA runs are:

| GEO Sample | SRA Run | Approximate SRA Size |
|---|---|---:|
| GSM8848584 | SRR32708522 | 1.38 GB |
| GSM8848585 | SRR32708521 | 1.03 GB |
| GSM8848586 | SRR32708520 | 185.69 MB |

The SRA files alone exceed the practical CellFlowX project budget before FASTQ conversion, temporary files, reference resources, and downstream analysis are considered.

## Decision

CellFlowX will begin from the author-provided processed Cell Ranger filtered feature-barcode matrices.

Expected input components include:

- `matrix.mtx.gz`
- `features.tsv.gz`
- `barcodes.tsv.gz`

These contain the gene-by-cell count matrices required for downstream single-cell transcriptomic analysis.

## Processed Data Budget

Approximate total processed input:

**68.3 MB**

This provides a substantial reduction in storage requirements while retaining the information required for the planned downstream analysis.

## Pipeline Boundary

CellFlowX begins at:

**Processed gene-cell count matrices**

and performs:

**Input Validation**
→ **AnnData Construction**
→ **QC**
→ **Filtering**
→ **Normalization**
→ **Feature Selection**
→ **Dimensionality Reduction**
→ **Clustering**
→ **Cell Annotation**
→ **Tumor Heterogeneity Analysis**
→ **Differential Analysis**
→ **Pathway Analysis**
→ **Reporting**

The following upstream stages are outside the CellFlowX execution boundary:

- SRA download
- FASTQ generation
- read-level QC
- barcode processing
- read alignment
- transcript quantification
- Cell Ranger execution

## Reference Genome Consequence

Because CellFlowX does not repeat alignment and quantification, it does not require downloading:

- GRCh38 FASTA
- gene annotation GTF
- STAR index
- Cell Ranger reference transcriptome

The original processed matrices remain the immutable starting point of the project.

## Why This Is Scientifically Appropriate

The principal scientific objectives of CellFlowX concern cell-level and population-level transcriptomic analysis rather than development of a read-processing pipeline.

Processed count matrices retain the information required for:

- cell QC
- gene QC
- normalization
- highly variable gene selection
- PCA
- neighborhood construction
- clustering
- UMAP
- marker identification
- cell-type annotation
- cell composition analysis
- transcriptional-state analysis
- differential analysis
- pathway analysis

## Engineering Consequences

### Advantages

- Fits local hardware constraints
- Keeps the project below the storage ceiling
- Reduces unnecessary compute
- Avoids duplicating upstream processing already performed by the study authors
- Allows engineering effort to focus on reproducible single-cell downstream workflows
- Makes CI/testing and local development more practical

### Trade-offs

CellFlowX cannot independently evaluate:

- raw sequencing quality
- alignment quality
- barcode calling
- Cell Ranger mapping performance
- upstream quantification parameters

These limitations will be documented rather than hidden.

## Reproducibility Strategy

Reproducibility will begin from the frozen processed matrices.

CellFlowX will record:

- GEO accessions
- sample identifiers
- input filenames
- file sizes
- checksums
- download provenance
- workflow parameters
- software versions

This defines a reproducible and auditable starting point without requiring raw sequencing reprocessing.

## Change Control

Raw FASTQ processing will not be added to CellFlowX unless the project scope and storage constraints are formally changed.

Any such change requires a new ADR.
