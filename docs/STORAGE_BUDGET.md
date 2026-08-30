# CellFlowX Storage Budget

## Storage Policy

CellFlowX is designed as a lightweight, local-first single-cell RNA-seq project.

The complete working project should remain approximately:

**≤ 1.5 GB**

This limit includes downloaded data, AnnData objects, analysis outputs, plots, tables, dashboard assets, and workflow files.

---

## Input Data Budget

CellFlowX uses processed 10x Genomics filtered feature-barcode matrices.

| GEO Sample | Approximate Download |
|---|---:|
| GSM8848584 | 32.1 MB |
| GSM8848585 | 31.4 MB |
| GSM8848586 | 4.8 MB |
| **Total** | **68.3 MB** |

Raw SRA/FASTQ files are excluded from the project.

---

## Estimated Project Budget

| Component | Target |
|---|---:|
| Downloaded processed matrices | ~70 MB |
| Extracted 10x matrices | ~100–200 MB |
| Working AnnData objects | ~200–500 MB |
| Tables and analysis results | ~50–100 MB |
| Figures | ~20–60 MB |
| Dashboard assets/data | ~30–100 MB |
| Code, tests and documentation | <30 MB |
| Temporary working allowance | ~200–300 MB |
| **Expected project footprint** | **~700 MB–1.2 GB** |
| **Hard ceiling** | **1.5 GB** |

These are planning estimates. Actual storage usage will be measured during execution.

---

## Reference Data

No additional reference genome resources are required.

CellFlowX will NOT download:

- GRCh38 FASTA
- GTF annotation
- STAR index
- Cell Ranger reference bundle

The selected input matrices have already undergone upstream alignment and gene quantification.

---

## AnnData Retention Policy

Large redundant `.h5ad` objects will not be permanently retained.

The preferred persistent objects are:

### 1. QC-filtered object

`data/processed/cellflowx_qc.h5ad`

Contains validated and filtered cells before downstream analysis.

### 2. Final annotated object

`data/processed/cellflowx_final.h5ad`

Contains the final analysis state required for downstream results and dashboard generation.

Temporary intermediate AnnData objects should be removed after validation when they are reproducible from the workflow.

---

## Raw Data Retention

Downloaded processed 10x matrices will be retained as the immutable project input.

They must never be manually modified.

Derived files must be written to separate processed or results directories.

---

## Temporary Files

Temporary files generated during:

- matrix extraction
- preprocessing
- Nextflow execution
- testing
- benchmarking

may be deleted after successful validation.

Nextflow `work/` will not be permanently retained after final reproducibility checks.

---

## GitHub Policy

Large biological data files must NOT be committed to GitHub.

Excluded files include:

- `.h5ad`
- `.h5`
- `.loom`
- `.rds`
- FASTQ
- BAM
- matrix archives
- generated Nextflow work directories

GitHub will contain:

- source code
- workflow definitions
- configuration
- manifests
- documentation
- tests
- lightweight result summaries
- selected lightweight figures where appropriate

---

## Storage Monitoring

Storage usage will be checked at major milestones using:

```bash
du -sh .
du -sh data/*
du -sh results/*