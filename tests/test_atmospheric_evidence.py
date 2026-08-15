import csv
from pathlib import Path


def test_atmospheric_evidence_is_source_graded():
    path = Path(__file__).resolve().parents[1] / "data" / "atmospheric_evidence.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows
    assert all(row["source_url"].startswith("https://") for row in rows)
    oxygen = [row for row in rows if row["species"] == "O2"]
    assert oxygen and all(row["status"] == "no evidence" for row in oxygen)
