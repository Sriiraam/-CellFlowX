from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
QC_CONFIG = ROOT / "config" / "qc_thresholds.json"


def test_qc_config_exists():
    assert QC_CONFIG.exists()


def test_qc_config_samples():
    config = json.loads(QC_CONFIG.read_text())

    expected = {
        "GSM8848584",
        "GSM8848585",
        "GSM8848586",
    }

    assert set(config) == expected


def test_qc_thresholds_valid():
    config = json.loads(QC_CONFIG.read_text())

    for sample, values in config.items():
        assert values["min_genes"] > 0
        assert values["max_genes"] > values["min_genes"]
        assert 0 < values["max_mt"] <= 100
