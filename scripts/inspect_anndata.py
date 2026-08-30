from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "data" / "processed" / "cellflowx_raw_merged.h5ad"


def main():
    print(f"Loading: {INPUT_FILE}")

    adata = sc.read_h5ad(INPUT_FILE)

    print("\n=== BASIC STRUCTURE ===")
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")
    print(f"Shape: {adata.shape}")

    assert adata.n_obs == 8660
    assert adata.n_vars == 18082

    print("\n=== MATRIX ===")
    print(f"Matrix type: {type(adata.X)}")
    print(f"Matrix dtype: {adata.X.dtype}")
    print(f"Sparse: {sparse.issparse(adata.X)}")

    if sparse.issparse(adata.X):
        values = adata.X.data
    else:
        values = np.asarray(adata.X).ravel()

    print(f"Minimum count: {values.min()}")
    print(f"Maximum count: {values.max()}")

    if np.any(values < 0):
        raise ValueError("Negative values detected in raw count matrix.")

    if not np.allclose(values, np.round(values)):
        raise ValueError("Non-integer values detected in raw count matrix.")

    print("Raw count validation: PASS")

    print("\n=== CELL METADATA ===")
    print(adata.obs.columns.tolist())
    print()
    print(adata.obs["geo_accession"].value_counts().sort_index())

    required_obs = {
        "geo_accession",
        "sample_id",
        "condition",
        "tissue",
    }

    missing_obs = required_obs - set(adata.obs.columns)

    if missing_obs:
        raise ValueError(f"Missing obs columns: {missing_obs}")

    print("\n=== GENE METADATA ===")
    print(adata.var.columns.tolist())

    print(f"Unique Ensembl IDs: {adata.var_names.is_unique}")

    if not adata.var_names.is_unique:
        raise ValueError("Ensembl IDs are not unique.")

    # Scanpy normally stores symbols here when gene_ids are var_names
    if "gene_symbols" in adata.var.columns:
        print(
            f"Gene symbols available: "
            f"{adata.var['gene_symbols'].notna().sum():,}"
        )
        print("\nExample genes:")
        print(adata.var[["gene_symbols"]].head())
    else:
        print("WARNING: gene_symbols column not found.")

    print("\n=== CELL BARCODE VALIDATION ===")
    print(f"Globally unique cell IDs: {adata.obs_names.is_unique}")

    if not adata.obs_names.is_unique:
        raise ValueError("Duplicated cell identifiers detected.")

    print("\n=== DATA LAYERS ===")
    print("Layers:", list(adata.layers.keys()))
    print("Raw attribute:", adata.raw)

    print("\n=================================")
    print("CellFlowX AnnData validation PASS")
    print("=================================")


if __name__ == "__main__":
    main()
