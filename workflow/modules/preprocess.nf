process PREPROCESS {

    tag "preprocess"

    publishDir params.processed_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_h5ad

    output:
    path "cellflowx_preprocessed.h5ad", emit: h5ad

    script:
    """
    python ${params.scripts_dir}/preprocess.py \
        --input ${input_h5ad} \
        --output cellflowx_preprocessed.h5ad
    """
}
