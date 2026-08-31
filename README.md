# 🧬 CellFlowX

[![CellFlowX CI](https://github.com/Sriiraam/CellFlowX/actions/workflows/ci.yml/badge.svg)](https://github.com/Sriiraam/CellFlowX/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.12-blue)
![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-23aa62)
![Scanpy](https://img.shields.io/badge/Scanpy-scRNA--seq-blueviolet)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED)
![License](https://img.shields.io/badge/License-MIT-green)

**CellFlowX** is a production-style, reproducible single-cell RNA-seq workflow for profiling **tumor heterogeneity in metastatic prostate cancer**.

It combines single-cell transcriptomic analysis with workflow engineering using **Scanpy, AnnData, Nextflow DSL2, Docker, SQLite, pytest, GitHub Actions, and Streamlit**.

---

## 🌐 Live Demo

Explore CellFlowX directly in your browser — no local installation required.

### **[Launch the CellFlowX Dashboard →](https://cellflowx.streamlit.app/)**

---

## 🧬 Workflow Architecture

![CellFlowX Workflow Architecture](docs/cellflowx_architecture.png)

The architecture summarizes the biological analysis, workflow orchestration, reproducibility, and deployment components of CellFlowX.

📊 **[View the complete CellFlowX Nextflow DAG →](docs/cellflowx_workflow_dag.html)**

---

## 🎯 Project Objective

CellFlowX investigates the following central biological question:

> **How do cellular composition and transcriptional states vary across metastatic prostate cancer samples, and which cell populations and molecular programs contribute to inter-tumor heterogeneity?**

The project focuses on:

- Single-cell quality control
- Sample-aware filtering
- Doublet detection
- Normalization and feature selection
- Dimensionality reduction
- Leiden clustering
- Cell-state annotation
- Inter-tumor heterogeneity
- CNV-like transcriptional evidence
- Tumor-state differential expression
- Functional enrichment
- Biological synthesis
- Workflow reproducibility
- Containerized execution
- Automated testing and CI/CD
- Interactive result exploration

---

## 📊 Dataset

CellFlowX uses publicly available single-cell RNA-seq data from:

**GEO:** `GSE292074`  
**BioProject:** `PRJNA1236646`  
**Organism:** Homo sapiens  
**Assay:** 10x Genomics single-cell RNA-seq  
**Study context:** Metastatic prostate cancer  
**Processing:** Cell Ranger v7 / GRCh38-2020-A

The implementation starts from **processed Cell Ranger filtered feature-barcode matrices** rather than downloading the substantially larger raw sequencing data.

### Samples

| GEO Sample | Sample ID | Initial Cells |
|---|---|---:|
| GSM8848584 | MH_07-042-M2 | 3,713 |
| GSM8848585 | MH_13-084-D12 | 3,988 |
| GSM8848586 | MH_13-084-D13 | 959 |
| **Total** | — | **8,660** |

Raw FASTQ/SRA files are not required for this implementation.

This design keeps the project computationally practical while retaining the complete processed single-cell expression matrices required for downstream analysis.

---

## 🔬 Workflow

CellFlowX is orchestrated as a **15-process Nextflow DSL2 workflow**.

```text
10x Filtered Matrices
        │
        ▼
AnnData Construction
        │
        ▼
QC Metrics
        │
        ├── QC Visualization
        │
        ▼
Sample-aware QC Evaluation
        │
        ▼
Doublet Detection
        │
        ▼
Final QC Filtering
        │
        ▼
Normalization + HVG Selection
        │
        ▼
PCA + Neighbors + Leiden + UMAP
        │
        ▼
Marker Discovery
        │
        ▼
Cell-state Annotation
        │
        ▼
Inter-tumor Heterogeneity
        │
        ▼
CNV-like Expression Evidence
        │
        ▼
Tumor-state Differential Expression
        │
        ▼
Functional Enrichment
        │
        ▼
Biological Synthesis
```

The workflow separates individual analytical stages into reusable Nextflow modules while maintaining a reproducible execution path from processed 10x matrices to biological interpretation.

📊 **[View the complete CellFlowX Nextflow DAG →](docs/cellflowx_workflow_dag.html)**

---

# 🔬 Biological Analysis

## 🧪 Quality Control

CellFlowX performs **sample-aware quality control** rather than applying one universal threshold across all three tumors.

QC evaluates:

- Total transcript counts
- Number of detected genes
- Mitochondrial transcript percentage
- Sample-specific filtering thresholds
- Scrublet doublet predictions
- Final cell retention

### Sample-specific QC thresholds

| Sample | Minimum Genes | Maximum Genes | Maximum MT % |
|---|---:|---:|---:|
| GSM8848584 | 500 | 9,500 | 20 |
| GSM8848585 | 500 | 9,000 | 8 |
| GSM8848586 | 450 | 6,500 | 8 |

These thresholds account for differences in the observed distributions between samples rather than assuming identical technical characteristics.

### Final QC

| Metric | Result |
|---|---:|
| Raw cells | 8,660 |
| Removed cells | 427 |
| Retained cells | 8,233 |
| Retention | 95.07% |
| Genes after preprocessing | 17,106 |
| Highly variable genes | 2,000 |

### Cells retained by sample

| Sample | QC-passed Cells |
|---|---:|
| GSM8848584 | 3,448 |
| GSM8848585 | 3,858 |
| GSM8848586 | 927 |
| **Total** | **8,233** |

Scrublet was executed independently for each sample.

A total of **7 predicted doublets** were conservatively excluded.

Ribosomal QC was not artificially introduced because canonical **RPL/RPS genes were absent from the supplied processed feature matrices**.

---

## 🧹 Preprocessing

After QC, CellFlowX performs:

1. Raw count preservation
2. Low-frequency gene filtering
3. Library-size normalization
4. Log transformation
5. Highly variable gene selection

Raw counts are preserved separately before normalization.

The final preprocessing stage contains:

- **8,233 cells**
- **17,106 genes**
- **2,000 highly variable genes**

The expression matrices remain sparse to avoid unnecessary memory expansion.

---

## 🗺️ Cell Atlas

CellFlowX constructs the single-cell representation using:

- Highly variable genes
- Principal component analysis
- Neighborhood graph construction
- Leiden community detection
- UMAP dimensionality reduction
- Cluster marker discovery

A total of **15 Leiden clusters** were identified.

### Major expression-state annotations

CellFlowX identified the following expression states:

- Prostate epithelial - AR high
- Prostate epithelial - luminal
- Neuroendocrine-like
- Steroidogenic-like
- Fibroblast
- Activated fibroblast
- Macrophage
- T/NK lymphocyte
- Endothelial
- Cycling
- Stress/hypoxia-like
- Epithelial-like - uncertain

These labels represent **expression-state annotations**.

They are not automatically interpreted as definitive malignant-cell identities.

---

## 🧬 Inter-Tumor Heterogeneity

CellFlowX compares cellular composition and transcriptional states across the three metastatic tumor samples.

Strong differences were observed between tumors.

| Sample | QC Cells | Dominant State | Proportion |
|---|---:|---|---:|
| GSM8848584 | 3,448 | Steroidogenic-like | 32.74% |
| GSM8848585 | 3,858 | Neuroendocrine-like | 43.47% |
| GSM8848586 | 927 | Activated fibroblast | 56.09% |

Additional prominent populations include:

- **GSM8848584:** AR-high prostate epithelial cells
- **GSM8848585:** Luminal prostate epithelial and cycling populations
- **GSM8848586:** Macrophage populations

The analysis therefore demonstrates substantial **inter-tumor cellular and transcriptional heterogeneity** across the three metastatic samples.

Strong sample-specific structure is also visible in the low-dimensional representation.

Because sample identity and biological state are partially confounded, these differences are interpreted carefully rather than automatically attributing them to a single biological mechanism.

---

## 🧬 CNV-like Malignancy Evidence

CellFlowX performs a chromosome-level expression-deviation analysis relative to reference cell populations.

The resulting score summarizes chromosome-scale transcriptional deviations and is used as **supporting evidence** when evaluating potentially malignant epithelial populations.

Strong CNV-like transcriptional evidence was observed particularly within:

- AR-high prostate epithelial populations
- Luminal prostate epithelial populations

Notable observations include:

| Population | Sample | CNV-like Positive Fraction |
|---|---|---:|
| Luminal prostate epithelial | GSM8848585 | 40.05% |
| AR-high prostate epithelial | GSM8848584 | 24.21% |

These measurements provide **supportive malignancy evidence only**.

They are:

- Not DNA-level CNV calls
- Not direct genomic measurements
- Not definitive proof of malignancy
- Not a replacement for validated DNA-based CNV analysis

Some stromal populations can also exhibit elevated scores because transcriptional variation and technical effects can influence chromosome-level expression patterns.

---

## 📈 Differential Expression

CellFlowX performs exploratory differential expression between major tumor-associated expression states.

### Comparisons

| Comparison | Significant Genes |
|---|---:|
| AR-high vs Neuroendocrine-like | 3,920 |
| Luminal vs Neuroendocrine-like | 3,112 |
| AR-high vs Steroidogenic-like | 4,255 |
| Neuroendocrine-like vs Steroidogenic-like | 4,429 |

Significance criteria:

```text
FDR < 0.05
|logFC| >= 1
```

These analyses are designed to identify genes distinguishing major transcriptional states.

However, they are interpreted as **exploratory cell-state differential expression** because several states are strongly associated with individual samples.

Cells are not treated as independent patient-level biological replicates.

---

## 🧠 Functional Enrichment

Significant transcriptional differences are further examined through functional enrichment analysis.

CellFlowX uses **GSEApy / Enrichr** with Gene Ontology Biological Process resources.

Representative biological programs include:

### Neuroendocrine-like

Enrichment includes programs associated with:

- Nervous-system development
- Chemical synaptic transmission

### Steroidogenic-like

Prominent enrichment includes:

- Cholesterol metabolism

### AR-high prostate epithelial

Comparisons involving AR-high populations highlight:

- Cholesterol-associated metabolic programs

### Luminal prostate epithelial

Prominent programs include:

- Regulation of cell migration

Functional enrichment is interpreted descriptively because only three biological samples are available and several expression states are sample-associated.

---

## 🧩 Biological Synthesis

The final biological synthesis integrates:

- QC results
- Cell-state annotation
- Sample composition
- Inter-tumor heterogeneity
- CNV-like transcriptional evidence
- Differential expression
- Functional enrichment

Together, these analyses demonstrate that the metastatic prostate cancer samples contain markedly different cellular compositions and transcriptional programs.

CellFlowX therefore provides an integrated view of **tumor-state diversity and inter-tumor heterogeneity** while explicitly retaining the limitations imposed by the small number of biological samples.

---

# ⚙️ Bioinformatics Engineering

## 🏗️ Engineering Architecture

CellFlowX combines biological analysis with production-style workflow engineering.

| Component | Technology |
|---|---|
| Workflow orchestration | Nextflow DSL2 |
| Single-cell analysis | Scanpy / AnnData |
| Programming | Python |
| Data representation | Sparse AnnData / H5AD |
| Containers | Docker |
| Database | SQLite |
| Testing | pytest |
| CI/CD | GitHub Actions |
| Dashboard | Streamlit |
| Version control | Git / GitHub |
| Provenance | Automated environment/run capture |
| Benchmarking | Nextflow trace/report summaries |

The architecture separates biological logic from workflow orchestration.

Python scripts implement individual analytical operations, while **Nextflow controls execution, dependencies, inputs, outputs, and portability**.

---

## 🔁 Nextflow DSL2

CellFlowX is organized into modular Nextflow DSL2 processes.

The production workflow contains **15 processes**:

```text
BUILD_ANNDATA
QC_METRICS
QC_PLOTS
EVALUATE_QC
DETECT_DOUBLETS
FINALIZE_QC
PREPROCESS
EMBEDDING_CLUSTERING
MARKER_DISCOVERY
ANNOTATE_CELLS
TUMOR_HETEROGENEITY
CNV_MALIGNANCY
TUMOR_STATE_DE
FUNCTIONAL_ENRICHMENT
BIOLOGICAL_SYNTHESIS
```

Individual modules are maintained under:

```text
workflow/modules/
```

This modular structure makes analytical stages easier to test, modify, resume, and reuse.

---

## 🚀 Execution Profiles

CellFlowX includes multiple Nextflow execution profiles.

Available profiles:

- `local`
- `docker`
- `slurm`
- `azure`
- `k8s-demo`

These profiles demonstrate how the workflow can move from local development toward containerized, HPC, cloud, and Kubernetes-style execution environments.

### Local execution

```bash
nextflow run main.nf -profile local
```

### Docker execution

```bash
nextflow run main.nf -profile docker
```

### Resume an interrupted workflow

```bash
nextflow run main.nf -profile docker -resume
```

Nextflow caching allows successfully completed processes to be reused when appropriate.

---

## 🐳 Docker

CellFlowX provides a containerized runtime environment to reduce dependency differences between systems.

The Docker environment includes the Python and system dependencies required by the analytical workflow.

The complete CellFlowX workflow was successfully validated through the Docker profile.

### Docker validation

All **15 processes completed successfully**:

```text
BUILD_ANNDATA              ✔
QC_METRICS                 ✔
QC_PLOTS                   ✔
EVALUATE_QC                ✔
DETECT_DOUBLETS            ✔
FINALIZE_QC                ✔
PREPROCESS                 ✔
EMBEDDING_CLUSTERING       ✔
MARKER_DISCOVERY           ✔
ANNOTATE_CELLS             ✔
TUMOR_HETEROGENEITY        ✔
CNV_MALIGNANCY             ✔
TUMOR_STATE_DE             ✔
FUNCTIONAL_ENRICHMENT      ✔
BIOLOGICAL_SYNTHESIS       ✔
```

The successful end-to-end Docker execution validates portability of the workflow beyond the original local Python environment.

---

## ⚡ Benchmarking

CellFlowX captures Nextflow execution metrics and converts them into compact benchmark summaries.

Generated benchmark artifacts include:

- `nextflow_report.html`
- `nextflow_timeline.html`
- `nextflow_trace.txt`
- `nextflow_dag.html`
- `process_benchmark_summary.csv`
- `project_benchmark_summary.csv`
- `BENCHMARKS.md`

### Project benchmark summary

| Metric | Result |
|---|---:|
| Tasks recorded | 15 |
| Processes recorded | 15 |
| Summed task runtime | 9.42 min |
| Peak task memory | 1.3 GB |
| Mean task CPU | 124.29% |

CPU utilization above 100% reflects multi-core process execution.

The **9.42-minute value represents summed task runtime**, not total pipeline wall-clock runtime.

Benchmarking provides a reproducible record of computational behavior in addition to biological results.

---

## 🗄️ SQLite Database

CellFlowX stores compact analytical summaries in a SQLite database.

### Database tables

```text
samples
celltype_composition
heterogeneity_summary
cnv_summary
de_summary
enrichment_summary
```

This provides a lightweight structured layer for dashboard summaries and analytical result retrieval.

Large cell-level expression matrices remain in **AnnData/H5AD** format.

This creates a hybrid architecture:

```text
H5AD
 │
 ├── Cell-level expression data
 ├── Embeddings
 ├── Cluster assignments
 └── Detailed metadata

SQLite
 │
 ├── Sample summaries
 ├── Cell-state composition
 ├── Heterogeneity summaries
 ├── CNV summaries
 ├── DE summaries
 └── Enrichment summaries
```

---

## 🖥️ Interactive Streamlit Dashboard

CellFlowX includes a multi-page Streamlit dashboard for interactive exploration of biological and engineering results.

### 🌐 Live application

**[Launch the CellFlowX Dashboard →](https://cellflowx.streamlit.app/)**

### Dashboard sections

1. Quality Control
2. Cell Atlas
3. Heterogeneity
4. CNV Evidence
5. Differential Expression
6. Pathway Enrichment
7. Biological Summary
8. Engineering & Performance

The dashboard combines:

- QC metrics
- UMAP exploration
- Cell-state composition
- Inter-tumor comparisons
- CNV-like evidence
- Differential expression
- Pathway enrichment
- Biological synthesis
- Workflow benchmarking

Lightweight deployment-specific data files are used so that large H5AD matrices do not need to be stored in the GitHub repository or deployed directly with the web application.

---

## ✅ Testing

CellFlowX includes automated repository validation using **pytest**.

Run local tests with:

```bash
pytest -q
```

Project-level validation can be executed using:

```bash
python scripts/validate_project.py
```

Validation covers repository components such as:

- Sample manifest structure
- Configuration
- SQLite database
- Processed biological outputs
- Repository integrity

Large processed biological datasets are intentionally excluded from GitHub and therefore remain part of **local validation rather than repository CI**.

---

## 🔄 CI/CD

CellFlowX uses **GitHub Actions** for automated repository checks.

The CI workflow automatically validates:

- Repository-safe pytest tests
- Python script compilation
- Nextflow configuration syntax

The workflow uses Python 3.12 and installs the repository dependencies before running validation.

The CI status is displayed directly at the top of this README.

---

## 🧾 Reproducibility & Provenance

CellFlowX records computational provenance for reproducibility.

Captured information includes:

- UTC execution timestamp
- Python version
- Operating system/platform
- Dependency versions
- Git commit
- Git branch
- Nextflow version
- Java version
- Docker information
- Dataset identifiers
- QC statistics
- Workflow configuration

Dependencies are frozen in:

```text
requirements.txt
```

QC thresholds are centralized in:

```text
config/qc_thresholds.json
```

Provenance outputs are stored separately from large biological result files so that reproducibility metadata can remain version controlled.

---

# 📁 Repository Structure

```text
CellFlowX/
│
├── main.nf
├── nextflow.config
├── requirements.txt
├── Dockerfile
│
├── config/
│   └── qc_thresholds.json
│
├── manifests/
│   └── samplesheet.csv
│
├── workflow/
│   ├── modules/
│   └── subworkflows/
│
├── scripts/
│
├── database/
│   ├── schema.sql
│   └── cellflowx.db
│
├── streamlit_app/
│   ├── app.py
│   ├── requirements.txt
│   ├── pages/
│   ├── utils/
│   └── data/
│
├── tests/
│
├── docs/
│   ├── cellflowx_architecture.png
│   └── cellflowx_workflow_dag.html
│
├── results/
│   ├── qc/
│   ├── annotation/
│   ├── heterogeneity/
│   ├── cnv/
│   ├── differential_expression/
│   ├── enrichment/
│   ├── synthesis/
│   ├── benchmarks/
│   └── provenance/
│
├── .github/
│   └── workflows/
│       └── ci.yml
│
├── LICENSE
├── CITATION.cff
└── README.md
```

---

# 🔄 Reproducible Analysis Strategy

CellFlowX follows several engineering principles intended to make the project easier to reproduce and maintain.

### 1. Workflow orchestration

Nextflow DSL2 manages analytical dependencies and execution order.

### 2. Modular analysis

Individual Python scripts implement specific biological operations.

### 3. Centralized configuration

QC thresholds and workflow parameters are maintained outside analytical scripts where appropriate.

### 4. Sparse single-cell matrices

Expression matrices remain sparse to minimize memory consumption.

### 5. Containerized execution

Docker provides a reproducible runtime environment.

### 6. Automated testing

pytest validates critical repository components.

### 7. Continuous integration

GitHub Actions checks repository integrity after changes.

### 8. Provenance capture

Software versions, Git state, configuration, and dataset identifiers are recorded.

### 9. Benchmarking

Nextflow execution metrics are converted into process and project summaries.

### 10. Lightweight deployment

The Streamlit deployment uses compact CSV/SQLite summaries instead of requiring the complete local H5AD dataset.

---

# ⚠️ Scientific Limitations

CellFlowX intentionally maintains conservative biological interpretation.

The dataset contains only **three metastatic tumor samples** and does not contain a normal/control cohort.

Therefore, CellFlowX does **not** make claims about:

- Tumor versus normal differences
- Primary versus metastatic disease
- Treatment response
- Survival
- Causality
- Patient-independent population effects

Cells are not considered independent biological replicates for patient-level inference.

The small number of biological samples limits population-level statistical conclusions.

Strong sample-specific structure is present in the dataset, meaning that some transcriptional states and differential-expression results are confounded with sample identity.

Automatic batch correction was therefore not used simply to remove sample-specific biological structure.

Expression-state annotations such as:

- Neuroendocrine-like
- AR-high prostate epithelial
- Luminal prostate epithelial

should not automatically be interpreted as definitive malignant identities.

CNV-like expression scores represent **chromosome-level transcriptional deviation** and are supporting evidence only.

They are not equivalent to:

- DNA sequencing
- Validated genomic CNV calling
- Cytogenetic measurements
- Definitive malignant-cell classification

Differential expression and functional enrichment are therefore interpreted primarily as **exploratory tumor-state analyses**.

---

# 📚 Data Availability

CellFlowX uses publicly available data from:

**GEO accession:** `GSE292074`  
**BioProject:** `PRJNA1236646`

The repository does not redistribute the complete biological dataset.

Large processed biological files and intermediate workflow outputs are intentionally excluded from Git version control.

The public repository instead contains:

- Workflow code
- Analysis scripts
- Configuration
- Sample manifest
- Tests
- Documentation
- Lightweight dashboard data
- Reproducibility metadata
- Benchmark summaries

---

# 🛠️ Technology Stack

| Area | Technology |
|---|---|
| Language | Python |
| Single-cell framework | Scanpy |
| Data structure | AnnData |
| Workflow engine | Nextflow DSL2 |
| Containers | Docker |
| Database | SQLite |
| Dashboard | Streamlit |
| Testing | pytest |
| CI/CD | GitHub Actions |
| Version control | Git / GitHub |
| Differential expression | Scanpy |
| Functional enrichment | GSEApy / Enrichr |
| Doublet detection | Scrublet |
| Visualization | Matplotlib / Plotly |
| Execution targets | Local / Docker / SLURM / Azure / Kubernetes demo |

---

# 🎯 Portfolio Focus

CellFlowX demonstrates both **single-cell bioinformatics** and **bioinformatics pipeline engineering**.

The project highlights experience with:

### Computational biology

- scRNA-seq analysis
- Single-cell QC
- Clustering
- Marker discovery
- Cell-state annotation
- Tumor heterogeneity
- Differential expression
- Functional enrichment
- Cancer transcriptomics

### Pipeline engineering

- Nextflow DSL2
- Modular workflow design
- Configuration management
- Containerization
- Reproducibility
- Provenance
- Benchmarking
- Testing
- CI/CD

### Data engineering and presentation

- AnnData/H5AD
- Sparse matrices
- SQLite
- Streamlit
- Interactive visualization
- Lightweight deployment architecture

---

# 👨‍💻 Author

**Sriram B**

B.Tech Biotechnology  
Bioinformatics / Computational Genomics

Project focus:

**Bioinformatics Pipeline Engineering · Single-Cell Genomics · Reproducible Computational Biology**

---

# 📄 License

CellFlowX is released under the **MIT License**.

See:

```text
LICENSE
```

for the complete license terms.

---

# 📖 Citation

Citation metadata are provided in:

```text
CITATION.cff
```

Software citation:

> **CellFlowX: Single-Cell Transcriptomic Profiling of Tumor Heterogeneity in Metastatic Prostate Cancer**

Version: **1.0.0**

---

# 🧬 Project Status

## **CellFlowX v1.0.0**

**Production-style portfolio release**

Current release includes:

- ✅ Complete biological workflow
- ✅ 15-process Nextflow DSL2 orchestration
- ✅ Sample-aware QC
- ✅ Doublet detection
- ✅ Single-cell clustering
- ✅ Cell-state annotation
- ✅ Inter-tumor heterogeneity analysis
- ✅ CNV-like transcriptional evidence
- ✅ Differential expression
- ✅ Functional enrichment
- ✅ Biological synthesis
- ✅ SQLite analytical summaries
- ✅ Interactive Streamlit dashboard
- ✅ Docker execution
- ✅ Local / Docker / SLURM / Azure / Kubernetes-demo profiles
- ✅ Automated testing
- ✅ GitHub Actions CI
- ✅ Benchmark reporting
- ✅ Provenance capture
- ✅ Workflow DAG
- ✅ Architecture documentation

### 🌐 Live Dashboard

**[Explore CellFlowX →](https://cellflowx.streamlit.app/)**

---

**CellFlowX — reproducible single-cell transcriptomics meets production-style bioinformatics engineering.**