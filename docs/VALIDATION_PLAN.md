
CellFlowX Validation Plan
Objective

Every major CellFlowX stage must produce technically valid outputs
before downstream execution continues.

Input Validation

Confirm:

all expected samples exist
required matrix files exist
matrix dimensions are valid
barcode count matches matrix columns
feature count matches matrix rows
sample identifiers are unique
gene identifiers are valid
AnnData Validation

Confirm:

cells are observations
genes are variables
sample metadata are present
count matrix is sparse where appropriate
raw counts are preserved
QC Validation

Inspect:

counts per cell
genes per cell
mitochondrial percentage
ribosomal percentage
doublet scores
cells retained per sample

Filtering decisions must be documented.

Preprocessing Validation

Confirm:

normalization completed successfully
log transformation is applied only where intended
highly variable genes are identified
raw counts required for later analysis remain available
Clustering Validation

Inspect:

PCA structure
neighbor graph
Leiden clusters
UMAP
sample distribution
cluster sizes
Annotation Validation

Cell-type labels must be supported by marker expression.

Ambiguous populations will remain conservatively labelled rather
than forced into unsupported cell types.

Differential Analysis Validation

Confirm:

comparison groups are explicitly defined
statistical unit is appropriate
raw/count data are used where required
limitations from sample number are reported
Output Validation

Required output files must exist and contain expected columns,
dimensions and metadata.

Automated Testing

pytest tests will cover critical reusable functions and input
validation logic.
