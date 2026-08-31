from pathlib import Path
import sqlite3

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "database" / "cellflowx.db"


def connect():
    return sqlite3.connect(DB)


def test_database_exists():
    assert DB.exists()


def test_database_integrity():
    with connect() as conn:
        result = conn.execute(
            "PRAGMA integrity_check;"
        ).fetchone()[0]

    assert result == "ok"


def test_required_tables():
    expected = {
        "samples",
        "celltype_composition",
        "heterogeneity_summary",
        "cnv_summary",
        "de_summary",
        "enrichment_summary",
    }

    with connect() as conn:
        tables = {
            row[0]
            for row in conn.execute(
                """
                SELECT name
                FROM sqlite_master
                WHERE type='table'
                """
            )
        }

    assert expected.issubset(tables)


def test_database_row_counts():
    expected = {
        "samples": 3,
        "celltype_composition": 36,
        "heterogeneity_summary": 3,
        "cnv_summary": 28,
        "de_summary": 4,
        "enrichment_summary": 8,
    }

    with connect() as conn:
        for table, count in expected.items():
            observed = conn.execute(
                f"SELECT COUNT(*) FROM {table}"
            ).fetchone()[0]

            assert observed == count
