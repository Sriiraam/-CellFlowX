process MARKER_DISCOVERY {

    tag "marker_discovery"

    publishDir params.annotation_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad

    output:
    path "cluster_markers.csv", emit: markers
    path "top_cluster_markers.csv", emit: top_markers
    path "canonical_marker_dotplot.png", emit: dotplot

    script:
    """
    python ${projectDir}/scripts/marker_discovery.py \
        --input "${input_h5ad}" \
        --markers "cluster_markers.csv" \
        --top-markers "top_cluster_markers.csv" \
        --dotplot "canonical_marker_dotplot.png"
    """
}
