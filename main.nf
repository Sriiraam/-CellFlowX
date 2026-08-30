nextflow.enable.dsl=2

include { BUILD_ANNDATA } from './workflow/modules/build_anndata'
include { QC_METRICS } from './workflow/modules/qc_metrics'
include { QC_PLOTS } from './workflow/modules/qc_plots'
include { EVALUATE_QC } from './workflow/modules/evaluate_qc'
include { DETECT_DOUBLETS } from './workflow/modules/detect_doublets'
include { FINALIZE_QC } from './workflow/modules/finalize_qc'
include { PREPROCESS } from './workflow/modules/preprocess'
include { EMBEDDING_CLUSTERING } from './workflow/modules/embedding_clustering'
include { MARKER_DISCOVERY } from './workflow/modules/marker_discovery'
include { ANNOTATE_CELLS } from './workflow/modules/annotate_cells'
include { TUMOR_HETEROGENEITY } from './workflow/modules/tumor_heterogeneity'
include { CNV_MALIGNANCY } from './workflow/modules/cnv_malignancy'
include { TUMOR_STATE_DE } from './workflow/modules/tumor_state_de'
include { FUNCTIONAL_ENRICHMENT } from './workflow/modules/functional_enrichment'
include { BIOLOGICAL_SYNTHESIS } from './workflow/modules/biological_synthesis'

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
        "${projectDir}/config/celltype_annotations.json",
        checkIfExists: true
    )

    ANNOTATE_CELLS(
        EMBEDDING_CLUSTERING.out.embedding,
        annotation_file
    )

    TUMOR_HETEROGENEITY(
        ANNOTATE_CELLS.out.annotated
    )

    CNV_MALIGNANCY(
        ANNOTATE_CELLS.out.annotated
    )

    TUMOR_STATE_DE(
        CNV_MALIGNANCY.out.h5ad
    )

    FUNCTIONAL_ENRICHMENT(
        TUMOR_STATE_DE.out.results
    )

    BIOLOGICAL_SYNTHESIS(
        TUMOR_HETEROGENEITY.out.summary
            .mix(TUMOR_HETEROGENEITY.out.percentages)
            .collect(),

        CNV_MALIGNANCY.out.summary
            .collect(),

        FUNCTIONAL_ENRICHMENT.out.results
    )
}
