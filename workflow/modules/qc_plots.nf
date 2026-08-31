process QC_PLOTS {

    tag "qc_plots"

    publishDir params.qc_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad

    output:
    path "qc_plots", emit: plots

    script:
    """
    mkdir -p qc_plots

    python ${params.scripts_dir}/plot_qc_distributions.py \
        --input ${input_h5ad} \
        --outdir qc_plots
    """
}
