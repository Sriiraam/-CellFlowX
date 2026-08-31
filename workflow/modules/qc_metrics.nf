process QC_METRICS {

    tag "qc_metrics"

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

    output:
    path "cellflowx_qc_metrics.h5ad", emit: h5ad
    path "qc_summary_by_sample.csv", emit: summary

    script:
    """
    python ${params.scripts_dir}/calculate_qc_metrics.py \
        --input ${input_h5ad} \
        --output cellflowx_qc_metrics.h5ad \
        --summary qc_summary_by_sample.csv
    """
}
