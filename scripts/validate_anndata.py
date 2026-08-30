from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_ROOT = PROJECT_ROOT / "data" / "raw" / "extracted"
OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "cellflowx_raw_merged.h5ad"
)


SAMPLES = {
    "GSM8848584": {
        "sample_id": "MH_07-042-M2",
        "expected_cells": 3713,
    },
    "GSM8848585": {
        "sample_id": "MH_13-084-D12",
        "expected_cells": 3988,
    },
    "GSM8848586": {
        "sample_id": "MH_13-084-D13",
        "expected_cells": 959,
    },
}

EXPECTED_GENES = 18082


def main():

    matrices = []
    obs_tables = []
    reference_var = None
    reference_genes = None

    for geo, metadata in SAMPLES.items():

        print(f"\nLoading {geo}...")

        adata = sc.read_10x_mtx(
            INPUT_ROOT / geo,
            var_names="gene_ids",
            make_unique=False,
            cache=False,
        )

        assert adata.n_obs == metadata["expected_cells"]
        assert adata.n_vars == EXPECTED_GENES
        assert adata.var_names.is_unique
        assert adata.X is not None
        assert sparse.issparse(adata.X)

        # Verify identical gene ordering across samples
        if reference_genes is None:
            reference_genes = adata.var_names.copy()
            reference_var = adata.var.copy()
        else:
            if not adata.var_names.equals(reference_genes):
                raise ValueError(
                    f"{geo}: gene order does not match reference sample"
                )

        # Globally unique cell IDs
        adata.obs_names = [
            f"{geo}_{barcode}"
            for barcode in adata.obs_names
        ]

        obs = pd.DataFrame(
            index=adata.obs_names
        )

        obs["geo_accession"] = geo
        obs["sample_id"] = metadata["sample_id"]
        obs["condition"] = "metastatic_prostate_cancer"
        obs["tissue"] = "tumor_tissue"

        matrices.append(adata.X.tocsr())
        obs_tables.append(obs)

        print(
            f"{geo}: "
            f"{adata.n_obs:,} cells × {adata.n_vars:,} genes"
        )

    print("\nCombining sparse matrices...")

    combined_X = sparse.vstack(
        matrices,
        format="csr",
    )

    combined_obs = pd.concat(
        obs_tables,
        axis=0,
    )

    merged = ad.AnnData(
        X=combined_X,
        obs=combined_obs,
        var=reference_var.copy(),
    )

    # Final validation
    assert merged.shape == (8660, 18082)
    assert merged.X is not None
    assert sparse.issparse(merged.X)
    assert merged.obs_names.is_unique
    assert merged.var_names.is_unique

    merged.uns["project"] = "CellFlowX"
    merged.uns["dataset"] = "GSE292074"
    merged.uns["bioproject"] = "PRJNA1236646"
    merged.uns["input_type"] = (
        "processed_10x_filtered_feature_bc_matrix"
    )

    print("\nFinal shape:", merged.shape)
    print("X type:", type(merged.X))
    print("Layers:", list(merged.layers.keys()))

    merged.write_h5ad(
        OUTPUT_FILE,
        compression="gzip",
    )

    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()