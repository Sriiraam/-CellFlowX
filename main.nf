nextflow.enable.dsl=2

include { BUILD_ANNDATA }    from './workflow/modules/build_anndata'
include { QC_METRICS }       from './workflow/modules/qc_metrics'
include { QC_PLOTS }         from './workflow/modules/qc_plots'
include { EVALUATE_QC }      from './workflow/modules/evaluate_qc'
include { DETECT_DOUBLETS }  from './workflow/modules/detect_doublets'
include { FINALIZE_QC }      from './workflow/modules/finalize_qc'
include { PREPROCESS }       from './workflow/modules/preprocess'



include { EMBEDDING_CLUSTERING } from './workflow/modules/embedding_clustering'


include { MARKER_DISCOVERY } from './workflow/modules/marker_discovery'


include { ANNOTATE_CELLS } from './workflow/modules/annotate_cells'

workflow {

    input_root = file(
        params.input_dir,
        checkIfExists: true
    )

    samplesheet = file(
        params.samplesheet,
        checkIfExists: true
    )

    thresholds = file(
        params.qc_thresholds,
        checkIfExists: true
    )

    BUILD_ANNDATA(
        input_root,
        samplesheet
    )

    QC_METRICS(
        BUILD_ANNDATA.out.h5ad
    )

    /*
     * Independent QC branches
     */
    QC_PLOTS(
        QC_METRICS.out.h5ad
    )

    EVALUATE_QC(
        QC_METRICS.out.h5ad,
        thresholds
    )

    DETECT_DOUBLETS(
        QC_METRICS.out.h5ad
    )

    FINALIZE_QC(
        DETECT_DOUBLETS.out.h5ad,
        thresholds
    )

    PREPROCESS(
        FINALIZE_QC.out.h5ad
    )

    EMBEDDING_CLUSTERING(
        PREPROCESS.out
    )

    MARKER_DISCOVERY(
        EMBEDDING_CLUSTERING.out.embedding
    )

    annotation_file = file(
        "${projectDir}/config/celltype_annotations.json"
    )

    ANNOTATE_CELLS(
        EMBEDDING_CLUSTERING.out.embedding,
        annotation_file
    )
}
