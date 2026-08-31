process BIOLOGICAL_SYNTHESIS {

    tag "biological_synthesis"

    publishDir params.synthesis_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path heterogeneity_files
    path cnv_files
    path enrichment_files

    output:
    path "phase12/*", emit: results

    script:
    """
    mkdir -p heterogeneity_input cnv_input enrichment_input phase12

    cp ${heterogeneity_files} heterogeneity_input/
    cp ${cnv_files} cnv_input/
    cp ${enrichment_files} enrichment_input/

    python ${params.scripts_dir}/biological_synthesis.py \
        --heterogeneity-dir heterogeneity_input \
        --cnv-dir cnv_input \
        --enrichment-dir enrichment_input \
        --output-dir phase12
    """
}
