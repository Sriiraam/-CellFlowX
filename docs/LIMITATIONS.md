
CellFlowX Limitations
Cohort Size

CellFlowX contains three metastatic prostate cancer scRNA-seq
samples.

This limits sample-level statistical power.

Control Group

The selected cohort does not contain a normal prostate control.

Therefore the project cannot directly infer tumor-versus-normal
transcriptional changes.

Raw Sequencing

CellFlowX begins from processed Cell Ranger count matrices.

Consequently, the project does not independently assess:

FASTQ quality
sequencing quality
alignment performance
barcode calling
Cell Ranger quantification performance
Biological Replication

Individual cells are not equivalent to independent patient-level
biological replicates.

Large cell numbers must not be interpreted as large biological
sample numbers.

Annotation

Cell-type annotation depends on observed markers and available
biological evidence.

Ambiguous clusters may remain broadly annotated.

Malignancy Assignment

Malignant-cell identification must be supported by available
expression evidence.

CellFlowX will not label cells malignant solely because they
originate from tumor tissue.

Generalization

Findings from these three samples cannot automatically be
generalized to all metastatic prostate cancer patients.

Interpretation

CellFlowX is primarily a reproducible exploratory and comparative
single-cell analysis project.

Causal or clinical claims are outside its scope.
