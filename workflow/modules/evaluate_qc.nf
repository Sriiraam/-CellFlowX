process EVALUATE_QC {

    tag "evaluate_qc"

    publishDir params.qc_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad
    path thresholds

    output:
    path "qc_threshold_evaluation.csv", emit: evaluation

    script:
    """
    python ${params.scripts_dir}/evaluate_qc_thresholds.py \
        --input "${input_h5ad}" \
        --thresholds "${thresholds}" \
        --output "qc_threshold_evaluation.csv"
    """
}
