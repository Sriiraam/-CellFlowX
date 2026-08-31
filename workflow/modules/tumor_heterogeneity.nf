process TUMOR_HETEROGENEITY {

    tag "tumor_heterogeneity"

    publishDir params.heterogeneity_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad

    output:
    path "celltype_counts_by_sample.csv", emit: counts
    path "celltype_percentages_by_sample.csv", emit: percentages
    path "heterogeneity_summary.csv", emit: summary
    path "celltype_composition_stacked.png", emit: stacked
    path "celltype_composition_heatmap.png", emit: heatmap

    script:
    """
    python ${params.scripts_dir}/tumor_heterogeneity.py \
        --input "${input_h5ad}" \
        --counts "celltype_counts_by_sample.csv" \
        --percentages "celltype_percentages_by_sample.csv" \
        --summary "heterogeneity_summary.csv" \
        --stacked-plot "celltype_composition_stacked.png" \
        --heatmap "celltype_composition_heatmap.png"
    """
}
