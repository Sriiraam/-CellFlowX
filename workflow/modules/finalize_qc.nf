process FINALIZE_QC {

    tag "finalize_qc"

    publishDir params.processed_dir,
        pattern: "*.h5ad",
        mode: 'copy',
        overwrite: true

    publishDir params.qc_results_dir,
        pattern: "*.csv",
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad
    path thresholds

    output:
    path "cellflowx_qc.h5ad", emit: h5ad
    path "cellflowx_qc_flagged.h5ad", emit: flagged
    path "final_qc_summary.csv", emit: summary
    path "qc_removal_reasons.csv", emit: removals
    path "doublet_candidates.csv", emit: doublets

    script:
    """
    python ${projectDir}/scripts/finalize_qc.py \
        --input ${input_h5ad} \
        --thresholds ${thresholds} \
        --output cellflowx_qc.h5ad \
        --flagged-output cellflowx_qc_flagged.h5ad \
        --summary final_qc_summary.csv \
        --removals qc_removal_reasons.csv \
        --doublets doublet_candidates.csv
    """
}
