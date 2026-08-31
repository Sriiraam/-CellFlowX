from pathlib import Path
import sqlite3
import json
import pandas as pd
import anndata as ad


ROOT = Path(__file__).resolve().parents[1]

checks = []


def check(name, condition):
    checks.append((name, bool(condition)))

    icon = "✓" if condition else "✗"
    print(f"{icon} {name}")


# --------------------------------------------------
# Manifest
# --------------------------------------------------

manifest_path = ROOT / "manifests" / "samplesheet.csv"

check(
    "Samplesheet exists",
    manifest_path.exists()
)

manifest = pd.read_csv(manifest_path)

check(
    "Three biological samples",
    len(manifest) == 3
)

check(
    "Sample IDs unique",
    manifest["sample_id"].is_unique
)

check(
    "Raw manifest cells = 8660",
    manifest["n_cells"].sum() == 8660
)


# --------------------------------------------------
# QC config
# --------------------------------------------------

qc_path = ROOT / "config" / "qc_thresholds.json"

check(
    "QC configuration exists",
    qc_path.exists()
)

qc = json.loads(qc_path.read_text())

check(
    "QC thresholds cover all samples",
    len(qc) == 3
)


# --------------------------------------------------
# SQLite
# --------------------------------------------------

db_path = ROOT / "database" / "cellflowx.db"

check(
    "SQLite database exists",
    db_path.exists()
)

with sqlite3.connect(db_path) as conn:

    integrity = conn.execute(
        "PRAGMA integrity_check"
    ).fetchone()[0]

    check(
        "SQLite integrity",
        integrity == "ok"
    )

    sample_count = conn.execute(
        "SELECT COUNT(*) FROM samples"
    ).fetchone()[0]

    check(
        "SQLite contains 3 samples",
        sample_count == 3
    )


# --------------------------------------------------
# Processed AnnData
# --------------------------------------------------

adata_path = (
    ROOT /
    "data" /
    "processed" /
    "cellflowx_annotated.h5ad"
)

check(
    "Annotated AnnData exists",
    adata_path.exists()
)

adata = ad.read_h5ad(
    adata_path,
    backed="r"
)

try:

    check(
        "QC-retained cells = 8233",
        adata.n_obs == 8233
    )

    check(
        "UMAP embedding present",
        "X_umap" in adata.obsm
    )

    check(
        "PCA embedding present",
        "X_pca" in adata.obsm
    )

    check(
        "Leiden clustering present",
        "leiden" in adata.obs.columns
    )

finally:
    adata.file.close()


# --------------------------------------------------
# Final result
# --------------------------------------------------

failed = [
    name
    for name, passed in checks
    if not passed
]

print()
print(
    f"Validation: "
    f"{len(checks) - len(failed)}/{len(checks)} passed"
)

if failed:

    print("\nFAILED:")

    for name in failed:
        print(f" - {name}")

    raise SystemExit(1)

print("\nCellFlowX project validation PASSED.")
