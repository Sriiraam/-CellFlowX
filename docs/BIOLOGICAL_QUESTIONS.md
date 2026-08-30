# CellFlowX Biological Questions

## Primary Biological Question

How do cellular composition and transcriptional states vary across metastatic prostate cancer samples, and which cell populations and molecular programs contribute to inter-tumor heterogeneity?

---

## Biological Objective

To characterize the single-cell landscape of metastatic prostate cancer and determine how cellular populations and transcriptional programs differ across the three tumor samples.

---

## Question 1 — Cellular Composition

**What major cell populations are present in the metastatic prostate cancer samples?**

Cell populations will be identified using:

- unsupervised clustering
- cluster-specific marker genes
- established lineage markers
- biological interpretation

Potential populations may include:

- malignant/epithelial cells
- T cells
- B cells
- monocytes/macrophages
- endothelial cells
- fibroblast/stromal cells

Final annotations will depend on observed marker expression.

---

## Question 2 — Inter-Tumor Heterogeneity

**How does cellular composition differ among GSM8848584, GSM8848585, and GSM8848586?**

We will compare:

- cell-type proportions
- cluster proportions
- sample-enriched populations
- sample-depleted populations
- overall cellular landscape

This represents the central comparative component of CellFlowX.

---

## Question 3 — Malignant-Cell Heterogeneity

**Do malignant/epithelial cells exhibit distinct transcriptional states across metastatic tumor samples?**

We will investigate:

- malignant-cell subclusters
- marker genes
- sample distribution
- transcriptional programs
- within-tumor and between-tumor heterogeneity

---

## Question 4 — Tumor Microenvironment

**How does the tumor microenvironment vary across metastatic prostate cancer samples?**

Where supported by the data, we will investigate differences in:

- immune populations
- stromal populations
- endothelial populations
- myeloid populations
- lymphoid populations

---

## Question 5 — Differential Transcriptional Programs

**Which genes and transcriptional programs distinguish major cell populations or tumor-associated states?**

Analyses may include:

- cluster marker identification
- cell-type-specific comparisons
- sample-associated transcriptional differences
- malignant-cell state comparisons

Statistical interpretation will account for the limited number of biological samples.

---

## Question 6 — Functional Pathways

**Which biological pathways characterize the major transcriptional states identified in the tumors?**

Potential analyses include:

- Gene Ontology Biological Process
- Reactome
- MSigDB Hallmark gene sets

Pathways will only be interpreted when supported by observed gene-expression results.

---

# Working Hypotheses

## Hypothesis 1

Metastatic prostate cancer samples will demonstrate substantial inter-tumor cellular heterogeneity.

## Hypothesis 2

The relative abundance of malignant, immune, and stromal populations will differ between tumor samples.

## Hypothesis 3

Malignant cells will contain multiple transcriptionally distinct states rather than forming one homogeneous population.

## Hypothesis 4

Distinct malignant-cell states will be associated with different biological pathway activities.

---

# Statistical Interpretation

The dataset contains only three metastatic tumor samples.

Therefore:

- individual cells will not automatically be treated as independent biological replicates
- sample-level conclusions will be interpreted cautiously
- pseudobulk/sample-aware approaches will be considered where statistically appropriate
- exploratory cell-state comparisons will be clearly distinguished from strong population-level inference

---

# Claims Outside Project Scope

CellFlowX will NOT claim:

- tumor versus healthy differences
- primary versus metastatic differences
- treatment response
- therapy resistance
- survival associations
- causal relationships
- population-wide prostate cancer biomarkers

unless the available metadata and analysis directly support such conclusions.

---

# Final Biological Story

CellFlowX aims to reconstruct the cellular landscape of metastatic prostate cancer at single-cell resolution and determine how tumor-cell states and tumor-microenvironment composition contribute to heterogeneity across metastatic tumor samples.