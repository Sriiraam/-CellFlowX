process CNV_MALIGNANCY {

    tag "cnv_malignancy"

    publishDir params.processed_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "cellflowx_cnv.h5ad"

    publishDir params.cnv_results_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "*.csv"

    publishDir params.cnv_results_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "*.png"

    input:
    path input_h5ad

    output:
    path "cellflowx_cnv.h5ad", emit: h5ad
    path "cnv_cell_scores.csv", emit: scores
    path "cnv_summary.csv", emit: summary
    path "cnv_umap.png", emit: umap
    path "cnv_chromosome_heatmap.png", emit: heatmap

    script:
    """
    python ${projectDir}/scripts/cnv_malignancy.py \
        --input "${input_h5ad}" \
        --output "cellflowx_cnv.h5ad" \
        --cell-scores "cnv_cell_scores.csv" \
        --summary "cnv_summary.csv" \
        --umap "cnv_umap.png" \
        --heatmap "cnv_chromosome_heatmap.png"
    """
}
