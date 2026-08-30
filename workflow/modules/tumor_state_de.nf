process TUMOR_STATE_DE {

    tag "tumor_state_de"

    publishDir params.de_results_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad

    output:
    path "phase10/*", emit: results

    script:
    """
    mkdir -p phase10

    python ${projectDir}/scripts/tumor_state_de.py \
        --input "${input_h5ad}" \
        --output-dir "phase10"
    """
}
