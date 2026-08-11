import gzip
import importlib.util
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).parents[1] / "bin" / "prepare_bedmethyl.py"
SPEC = importlib.util.spec_from_file_location("prepare_bedmethyl", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def bedmethyl_row(start, strand=".", modified=5, chrom="chrX"):
    fields = [
        chrom, str(start), str(start + 1), "m", "10", strand,
        str(start), str(start + 1), "255,0,0", "10", "50.0",
        str(modified), "0", str(10 - modified), "0", "0", "0", "0",
    ]
    return "\t".join(fields) + "\n"


def write_bedmethyl(path, rows):
    with gzip.open(path, "wt") as handle:
        handle.writelines(rows)


def test_same_coordinate_strands_are_combined_for_dss(tmp_path):
    source = tmp_path / "sample.bedmethyl.gz"
    write_bedmethyl(source, [
        bedmethyl_row(73820656, "-", 4),
        bedmethyl_row(73820656, "+", 6),
        bedmethyl_row(73820657, ".", 5),
    ])

    records = MODULE.coalesced(MODULE.bedmethyl_records(source, "m", None))

    assert [(row[1], row[4], row[5], row[6]) for row in records] == [
        (73820656, ".", 20, 10),
        (73820657, ".", 10, 5),
    ]


def test_coordinate_regression_is_rejected(tmp_path):
    source = tmp_path / "sample.bedmethyl.gz"
    write_bedmethyl(source, [bedmethyl_row(20), bedmethyl_row(10)])

    with pytest.raises(SystemExit, match="not coordinate sorted"):
        list(MODULE.bedmethyl_records(source, "m", None))


def test_region_ignores_order_of_unselected_contigs(tmp_path):
    source = tmp_path / "sample.bedmethyl.gz"
    write_bedmethyl(source, [
        bedmethyl_row(200, chrom="chr2"),
        bedmethyl_row(100, chrom="chr1"),
        bedmethyl_row(73820656),
        bedmethyl_row(73820657),
    ])
    region = MODULE.parse_region("chrX:73800000-73900000")

    records = list(MODULE.bedmethyl_records(source, "m", region))

    assert [(row[0], row[1]) for row in records] == [
        ("chrX", 73820656),
        ("chrX", 73820657),
    ]
