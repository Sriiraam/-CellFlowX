# CellFlowX Data Provenance

## Study

- GEO: GSE292074
- BioProject: PRJNA1236646
- Organism: Homo sapiens
- Disease: Metastatic prostate cancer
- Assay: scRNA-seq
- Input: Processed 10x Genomics filtered feature-barcode matrices

## Selected Samples

| Sample | GEO | SRA | Cells | Genes |
|---|---|---|---:|---:|
| MH_07-042-M2 | GSM8848584 | SRR32708522 | 3713 | 18082 |
| MH_13-084-D12 | GSM8848585 | SRR32708521 | 3988 | 18082 |
| MH_13-084-D13 | GSM8848586 | SRR32708520 | 959 | 18082 |

Total cells: **8660**

## Downloaded Files

| GEO | File | Approx Size | SHA256 |
|---|---|---:|---|
| GSM8848584 | GSM8848584_MH_07-042-M2_filtered_feature_bc_matrix.tar.gz | 33 MB | 5864f341bef961935b20f131b153567e38b64947a664cfc4c950c5c2e36ff3fc |
| GSM8848585 | GSM8848585_MH_13-084-D12_filtered_feature_bc_matrix.tar.gz | 32 MB | 8cab6e275dd62026600c3a5c7133bf829e756543388d339ba89447eb5e901ac3 |
| GSM8848586 | GSM8848586_MH_13-084-D13_filtered_feature_bc_matrix.tar.gz | 4.8 MB | 0010619840944bf3b3cac4261bdd0e88c2bd4029f3a312b1b73ecdd400efa0d6 |

## Matrix Validation

All three samples contain:

- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz

Validation results:

- matrix rows match feature count
- matrix columns match barcode count
- all barcodes are unique
- all Ensembl IDs are unique
- all samples share the same 18,082-gene feature space

Gene symbols contain 3 duplicated names.

Therefore:

**Ensembl IDs will be used as the primary unique gene identifier.**

Gene symbols will be retained as annotation.

## Storage

Downloaded archives: ~69 MB

Extracted matrices: ~69 MB

Current biological-data footprint: ~138 MB
