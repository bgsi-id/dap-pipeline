#!/usr/bin/env python3
"""Write the cohort summary report from a resolved DAP input manifest.

usage: summarize_cohort.py <resolved-inputs.json> <report.tsv> <report.json>

Only the aggregate `summary` block of the manifest is used. Sample identifiers
in the manifest are deliberately ignored: this report is aggregate-only.
"""
import csv
import json
import sys

STATS = ("n", "missing", "mean", "q1", "median", "q3")
_FORMULA_PREFIXES = ("=", "+", "-", "@")


def _safe_cell(value):
    if isinstance(value, str) and value.startswith(_FORMULA_PREFIXES):
        try:
            float(value)
        except ValueError:
            return f"'{value}"
    return value


def report_rows(summary):
    min_cell = summary.get("min_cell", 5)
    yield ("overview", "samples", "total", summary.get("sample_count"))
    for field in summary.get("fields") or []:
        label = field.get("label") or field.get("key")
        for stat in STATS:
            if stat in field:
                yield ("demographics", label, stat, field[stat])
        for category, count in (field.get("counts") or {}).items():
            yield ("demographics", label, category, count)
    for role, count in (summary.get("data_availability") or {}).items():
        yield ("data_availability", role, "samples", count)
    yield ("notes", "small_counts", f"values below {min_cell} are shown as <{min_cell}", "")


def main(argv):
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 64
    manifest_path, tsv_path, json_path = argv[1:]
    with open(manifest_path) as handle:
        summary = json.load(handle).get("summary")
    if not isinstance(summary, dict):
        print(
            "The resolved input manifest has no summary; the data service did "
            "not return one for this pack.",
            file=sys.stderr,
        )
        return 65
    with open(tsv_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["section", "field", "item", "value"])
        for row in report_rows(summary):
            writer.writerow([_safe_cell(cell) for cell in row])
    with open(json_path, "w") as handle:
        json.dump({"summary": summary}, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
