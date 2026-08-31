from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RESULTS = PROJECT_ROOT / "streamlit_app" / "data" / "results"
DATA = PROJECT_ROOT / "data"

QC_DIR = RESULTS / "qc"
EMBEDDING_DIR = RESULTS / "embedding"
ANNOTATION_DIR = RESULTS / "annotation"
HETEROGENEITY_DIR = RESULTS / "heterogeneity"
CNV_DIR = RESULTS / "cnv"
DE_DIR = RESULTS / "differential_expression" / "phase10"
ENRICHMENT_DIR = RESULTS / "enrichment" / "phase11"
SYNTHESIS_DIR = RESULTS / "synthesis" / "phase12"

ANNOTATED_H5AD = DATA / "processed" / "cellflowx_annotated.h5ad"
CNV_H5AD = DATA / "processed" / "cellflowx_cnv.h5ad"
