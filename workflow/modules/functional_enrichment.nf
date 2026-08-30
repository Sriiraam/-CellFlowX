process FUNCTIONAL_ENRICHMENT {

    tag "functional_enrichment"

    publishDir params.enrichment_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path phase10_files

    output:
    path "phase11/*", emit: results

    script:
    """
    mkdir -p phase10_input phase11

    cp ${phase10_files} phase10_input/

    python ${projectDir}/scripts/functional_enrichment.py \
        --de-dir phase10_input \
        --output-dir phase11
    """
}
