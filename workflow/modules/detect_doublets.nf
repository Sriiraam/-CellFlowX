process DETECT_DOUBLETS {

    tag "detect_doublets"

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
    path "cellflowx_qc_doublets.h5ad", emit: h5ad
    path "doublet_summary.csv", emit: summary

    script:
    """
    python ${params.scripts_dir}/detect_doublets.py \
        --input ${input_h5ad} \
        --output cellflowx_qc_doublets.h5ad \
        --summary doublet_summary.csv
    """
}
