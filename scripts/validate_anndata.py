import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-root",
        required=True,
    )

    parser.add_argument(
        "--samplesheet",
        required=True,
    )

    parser.add_argument(
        "--output",
        required=True,
    )

    return parser.parse_args()


def main():
    args = parse_args()

    input_root = Path(args.input_root)
    samplesheet = pd.read_csv(args.samplesheet)
    output = Path(args.output)

    matrices = []
    obs_tables = []
    reference_var = None

    for _, row in samplesheet.iterrows():

        geo = row["geo_accession"]

        print(f"\nLoading {geo}...")

        sample_dir = input_root / geo

        adata = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_ids",
            make_unique=False,
        )

        print(
            f"{geo}: "
            f"{adata.n_obs:,} cells × "
            f"{adata.n_vars:,} genes"
        )

        if reference_var is None:
            reference_var = adata.var.copy()

        else:
            if not adata.var_names.equals(
                reference_var.index
            ):
                raise ValueError(
                    f"Gene ordering differs for {geo}"
                )

        adata.obs_names = [
            f"{geo}_{barcode}"
            for barcode in adata.obs_names
        ]

        obs = pd.DataFrame(
            index=adata.obs_names
        )

        obs["geo_accession"] = geo
        obs["sample_id"] = row["sample_id"]
        obs["condition"] = row["condition"]
        obs["tissue"] = row["tissue"]

        matrices.append(
            adata.X.tocsr()
        )

        obs_tables.append(obs)

    print("\nCombining sparse matrices...")

    combined_X = sparse.vstack(
        matrices,
        format="csr",
    )

    combined_obs = pd.concat(
        obs_tables
    )

    merged = ad.AnnData(
        X=combined_X,
        obs=combined_obs,
        var=reference_var.copy(),
    )

    if merged.X is None:
        raise ValueError(
            "Merged X is missing"
        )

    if not sparse.issparse(
        merged.X
    ):
        raise ValueError(
            "Merged matrix is not sparse"
        )

    if not merged.obs_names.is_unique:
        raise ValueError(
            "Cell IDs are not unique"
        )

    if not merged.var_names.is_unique:
        raise ValueError(
            "Gene IDs are not unique"
        )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    merged.write_h5ad(
        output,
        compression="gzip",
    )

    print(
        f"\nFinal shape: {merged.shape}"
    )

    print(
        f"X type: {type(merged.X)}"
    )

    print(
        f"Layers: "
        f"{list(merged.layers.keys())}"
    )

    print(
        f"Saved: {output.resolve()}"
    )


if __name__ == "__main__":
    main()
