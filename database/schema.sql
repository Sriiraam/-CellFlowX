DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS celltype_composition;
DROP TABLE IF EXISTS heterogeneity_summary;
DROP TABLE IF EXISTS cnv_summary;
DROP TABLE IF EXISTS de_summary;
DROP TABLE IF EXISTS enrichment_summary;

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    geo_accession TEXT,
    condition TEXT,
    tissue TEXT,
    organism TEXT,
    assay TEXT,
    n_cells INTEGER,
    n_genes INTEGER
);

CREATE TABLE celltype_composition (
    sample_id TEXT,
    cell_state TEXT,
    percentage REAL,
    PRIMARY KEY (sample_id, cell_state)
);

CREATE TABLE heterogeneity_summary (
    sample_id TEXT PRIMARY KEY,
    total_cells INTEGER,
    detected_states INTEGER,
    dominant_state TEXT,
    dominant_state_pct REAL
);

CREATE TABLE cnv_summary (
    sample_id TEXT,
    cell_state TEXT,
    cnv_high_pct REAL,
    median_cnv_score REAL,
    PRIMARY KEY (sample_id, cell_state)
);

CREATE TABLE de_summary (
    comparison TEXT PRIMARY KEY,
    state_a TEXT,
    state_b TEXT,
    cells_a INTEGER,
    cells_b INTEGER,
    significant_genes INTEGER,
    up_in_a INTEGER,
    up_in_b INTEGER
);

CREATE TABLE enrichment_summary (
    gene_set TEXT PRIMARY KEY,
    input_genes INTEGER,
    significant_pathways INTEGER,
    top_pathway TEXT,
    top_adjusted_p REAL
);
