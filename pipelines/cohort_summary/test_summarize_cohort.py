import csv
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent / "bin" / "summarize_cohort.py"

SUMMARY = {
    "min_cell": 5,
    "sample_count": 3469,
    "fields": [
        {"key": "age_at_recruitment", "label": "Age at recruitment", "type": "number",
         "n": 3469, "missing": 0, "mean": 33.4, "q1": 24.0, "median": 32.0, "q3": 42.0},
        {"key": "sex", "label": "Sex", "type": "string", "n": 3469, "missing": 0,
         "counts": {"female": 2400, "male": 1069}},
        {"key": "patient_categ", "label": "Patient category", "type": "string",
         "n": 2527, "missing": 942, "counts": {"a": 2500, "b": "<5", "c": "<5"}},
    ],
    "data_availability": {"vcf": 3334, "cnv": "<5"},
}


def run(tmp_path, manifest):
    source = tmp_path / "resolved-inputs.json"
    source.write_text(json.dumps(manifest))
    done = subprocess.run(
        [sys.executable, str(SCRIPT), str(source), str(tmp_path / "out.tsv"), str(tmp_path / "out.json")],
        capture_output=True,
        text=True,
    )
    return done


def rows(tmp_path):
    with open(tmp_path / "out.tsv", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def value(table, section, field, item):
    return next(r["value"] for r in table if (r["section"], r["field"], r["item"]) == (section, field, item))


def test_report_lists_overview_demographics_and_data_availability(tmp_path):
    done = run(tmp_path, {"summary": SUMMARY})

    assert done.returncode == 0, done.stderr
    table = rows(tmp_path)
    assert value(table, "overview", "samples", "total") == "3469"
    assert value(table, "demographics", "Sex", "female") == "2400"
    assert value(table, "demographics", "Age at recruitment", "median") == "32.0"
    assert value(table, "demographics", "Age at recruitment", "n") == "3469"
    assert value(table, "data_availability", "vcf", "samples") == "3334"


def test_suppressed_counts_stay_suppressed(tmp_path):
    run(tmp_path, {"summary": SUMMARY})

    table = rows(tmp_path)
    assert value(table, "demographics", "Patient category", "b") == "<5"
    assert value(table, "data_availability", "cnv", "samples") == "<5"
    assert value(table, "demographics", "Patient category", "missing") == "942"


def test_json_report_carries_the_same_summary(tmp_path):
    run(tmp_path, {"summary": SUMMARY})

    report = json.loads((tmp_path / "out.json").read_text())
    assert report["summary"] == SUMMARY


def test_sample_identifiers_never_reach_the_report(tmp_path):
    run(tmp_path, {"summary": SUMMARY, "samples": [{"id": "SECRET-SAMPLE-ID", "assets": []}]})

    for name in ("out.tsv", "out.json"):
        assert "SECRET-SAMPLE-ID" not in (tmp_path / name).read_text()


def test_a_manifest_without_a_summary_fails_loudly(tmp_path):
    done = run(tmp_path, {"samples": []})

    assert done.returncode != 0
    assert "no summary" in done.stderr.lower()
    assert not (tmp_path / "out.tsv").exists()
