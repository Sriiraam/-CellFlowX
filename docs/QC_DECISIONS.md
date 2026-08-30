# CellFlowX QC Decisions

Input cells: 8,660
Retained cells: 8,233
Removed cells: 427
Retention: 95.07%

Sample retention:
- GSM8848584: 3,448 / 3,713
- GSM8848585: 3,858 / 3,988
- GSM8848586: 927 / 959

QC metrics:
- total counts
- detected genes
- mitochondrial percentage
- doublet score

Canonical RPL/RPS genes were absent from the supplied processed feature matrix,
therefore ribosomal percentage was not used for filtering.

Sample-aware thresholds:
- GSM8848584: min genes 500, max genes 9500, max MT 20%
- GSM8848585: min genes 500, max genes 9000, max MT 8%
- GSM8848586: min genes 450, max genes 6500, max MT 8%

Doublet detection:
Scrublet was run independently per sample.
Predicted doublets:
- GSM8848584: 0
- GSM8848585: 4
- GSM8848586: 3

Final QC object:
data/processed/cellflowx_qc.h5ad

