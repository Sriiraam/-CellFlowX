from pathlib import Path
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]

ANNOTATED = (
    ROOT /
    "data" /
    "processed" /
    "cellflowx_annotated.h5ad"
)


def test_annotated_file_exists():
    assert ANNOTATED.exists()


def test_annotated_shape():
    adata = ad.read_h5ad(
        ANNOTATED,
        backed="r"
    )

    try:
        assert adata.n_obs == 8233
        assert adata.n_vars > 17000
    finally:
        adata.file.close()


def test_embedding_present():
    adata = ad.read_h5ad(
        ANNOTATED,
        backed="r"
    )

    try:
        assert "X_umap" in adata.obsm
        assert "X_pca" in adata.obsm
        assert "leiden" in adata.obs.columns
    finally:
        adata.file.close()


def test_annotation_present():
    adata = ad.read_h5ad(
        ANNOTATED,
        backed="r"
    )

    candidates = {
        "cell_type",
        "celltype",
        "cell_state",
        "annotation",
        "cell_type_annotation",
    }

    try:
        assert any(
            col in adata.obs.columns
            for col in candidates
        )
    finally:
        adata.file.close()
