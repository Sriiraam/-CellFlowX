process ANNOTATE_CELLS {

    tag "annotate_cells"

    publishDir params.processed_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "cellflowx_annotated.h5ad"

    publishDir params.annotation_results_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "*.csv"

    publishDir params.annotation_results_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "*.png"

    input:
    path input_h5ad
    path annotations

    output:
    path "cellflowx_annotated.h5ad", emit: annotated
    path "annotated_umap.png", emit: umap
    path "cluster_annotations.csv", emit: cluster_table
    path "sample_celltype_composition.csv", emit: sample_table

    script:
    """
    python ${projectDir}/scripts/annotate_cells.py \
        --input "${input_h5ad}" \
        --annotations "${annotations}" \
        --output "cellflowx_annotated.h5ad" \
        --umap "annotated_umap.png" \
        --cluster-table "cluster_annotations.csv" \
        --sample-table "sample_celltype_composition.csv"
    """
}
