# CellFlowX Dataset Decision

## Decision Status

**Status: FROZEN**

The dataset selected for CellFlowX is:

- GEO Series: GSE292074
- BioProject: PRJNA1236646
- Organism: Homo sapiens
- Disease: Metastatic prostate cancer
- Assay: Single-cell RNA sequencing
- Technology: 10x Genomics
- Source: Tumor tissue

## Selected Samples

| Sample ID | GEO Accession | Biological Source |
|---|---|---|
| MH_07-042-M2 | GSM8848584 | Metastatic prostate cancer tumor tissue |
| MH_13-084-D12 | GSM8848585 | Metastatic prostate cancer tumor tissue |
| MH_13-084-D13 | GSM8848586 | Metastatic prostate cancer tumor tissue |

All three available scRNA-seq samples from the selected study are included.

## Selected Input

CellFlowX will use the author-provided processed Cell Ranger filtered feature-barcode matrices.

Expected 10x matrix components include:

- matrix.mtx.gz
- features.tsv.gz
- barcodes.tsv.gz

These matrices contain gene-by-cell count data suitable for downstream analysis using Scanpy and AnnData.

## Download Budget

Approximate processed input:

| Sample | Approximate Size |
|---|---:|
| GSM8848584 | 32.1 MB |
| GSM8848585 | 31.4 MB |
| GSM8848586 | 4.8 MB |
| **Total** | **68.3 MB** |

The complete CellFlowX project has a target storage ceiling of approximately **1.5 GB**.

## Why Processed Counts Were Selected

Raw sequencing data are available through the Sequence Read Archive (SRA), but the raw runs are substantially larger than the processed matrices.

Using raw sequencing data would additionally require:

- FASTQ conversion
- 10x Genomics read processing
- reference genome resources
- Cell Ranger
- substantially more temporary storage
- substantially higher memory and compute requirements

This would not provide enough additional value for the primary objective of CellFlowX, which is downstream single-cell transcriptomic analysis and reproducible workflow engineering.

Therefore, processed count matrices are the appropriate project input.

## Raw Data Policy

The following SRA runs will NOT be downloaded for the CellFlowX workflow:

- SRR32708522 — GSM8848584
- SRR32708521 — GSM8848585
- SRR32708520 — GSM8848586

Raw FASTQ files are outside the frozen CellFlowX execution scope.

## Reference Genome Decision

No genome FASTA, GTF, STAR index, or Cell Ranger reference bundle will be downloaded.

The author-provided matrices have already undergone upstream read processing and gene quantification.

CellFlowX begins from the gene-count matrix stage.

## Biological Rationale

This dataset provides true human single-cell transcriptomic data from metastatic prostate cancer.

It enables CellFlowX to investigate:

- cellular composition
- tumor microenvironment structure
- inter-tumor heterogeneity
- malignant-cell transcriptional states
- sample-specific cell populations
- marker genes
- differential transcriptional programs
- pathway-level biological differences

## Important Limitation

The dataset does not provide a normal prostate control group within the selected scRNA-seq cohort.

Therefore CellFlowX will NOT make claims about:

- tumor versus healthy prostate
- cancer-specific changes relative to normal tissue
- treatment response

The primary comparison is heterogeneity among metastatic tumor samples.

## Decision

**GSE292074 is frozen as the CellFlowX dataset.**

The project will use all three selected processed scRNA-seq samples and will not arbitrarily downsample the study by selecting individual sequencing runs.

Any future change to this dataset decision must be explicitly documented through a new Architecture Decision Record (ADR).