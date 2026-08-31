from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "manifests" / "samplesheet.csv"


def test_manifest_exists():
    assert MANIFEST.exists()


def test_manifest_schema():
    df = pd.read_csv(MANIFEST)

    required = {
        "sample_id",
        "geo_accession",
        "condition",
        "tissue",
        "organism",
        "assay",
        "input_type",
        "n_cells",
        "n_genes",
    }

    assert required.issubset(df.columns)


def test_manifest_samples():
    df = pd.read_csv(MANIFEST)

    assert len(df) == 3
    assert df["sample_id"].is_unique
    assert df["geo_accession"].is_unique
    assert df["n_cells"].sum() == 8660
