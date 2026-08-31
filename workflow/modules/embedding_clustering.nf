process EMBEDDING_CLUSTERING {

    tag "embedding_clustering"

    publishDir params.processed_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "cellflowx_embedding.h5ad"

    publishDir params.embedding_results_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "*.csv"

    publishDir params.embedding_results_dir,
        mode: 'copy',
        overwrite: true,
        pattern: "*.png"

    input:
    path input_h5ad

    output:
    path "cellflowx_embedding.h5ad",
        emit: embedding

    path "pca_variance.csv",
        emit: pca_variance

    path "cluster_summary.csv",
        emit: cluster_summary

    path "sample_cluster_composition.csv",
        emit: sample_cluster

    path "pca_explained_variance.png",
        emit: pca_plot

    path "umap_by_sample.png",
        emit: umap_sample

    path "umap_by_cluster.png",
        emit: umap_cluster

    script:
    """
    python ${params.scripts_dir}/embedding_clustering.py \
        --input "${input_h5ad}" \
        --output "cellflowx_embedding.h5ad" \
        --pca-variance "pca_variance.csv" \
        --cluster-summary "cluster_summary.csv" \
        --sample-cluster-summary "sample_cluster_composition.csv" \
        --pca-plot "pca_explained_variance.png" \
        --umap-sample "umap_by_sample.png" \
        --umap-cluster "umap_by_cluster.png" \
        --n-pcs 30 \
        --n-neighbors 15 \
        --resolution 0.5
    """
}
