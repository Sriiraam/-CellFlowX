process BUILD_ANNDATA {

    tag "build_anndata"

    publishDir params.processed_dir,
        mode: 'copy',
        overwrite: true

    input:
    path input_root
    path samplesheet

    output:
    path "cellflowx_raw_merged.h5ad", emit: h5ad

    script:
    """
    python ${projectDir}/scripts/validate_anndata.py \
        --input-root ${input_root} \
        --samplesheet ${samplesheet} \
        --output cellflowx_raw_merged.h5ad
    """
}
