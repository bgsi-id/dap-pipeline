import base64
import gzip
import importlib.util
import sys
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).parents[1] / "bin" / "prepare_bedmethyl.py"
SPEC = importlib.util.spec_from_file_location("prepare_bedmethyl", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def bedmethyl_row(start, strand=".", modified=5, chrom="chrX", coverage=10):
    fields = [
        chrom, str(start), str(start + 1), "m", str(coverage), strand,
        str(start), str(start + 1), "255,0,0", str(coverage), f"{(modified/coverage)*100:.1f}",
        str(modified), "0", str(coverage - modified), "0", "0", "0", "0",
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


def test_methylong_pre_combined_strands_output(tmp_path):
    """
    Simulates nf-core/methylong output where strands are pre-combined (--combine-strands).
    Records have strand '.' and unique coordinates.
    """
    source = tmp_path / "sample.methylong.bedmethyl.gz"
    write_bedmethyl(source, [
        bedmethyl_row(100, ".", modified=3, chrom="chr1", coverage=20),
        bedmethyl_row(200, ".", modified=18, chrom="chr1", coverage=20),
    ])

    records = list(MODULE.coalesced(MODULE.bedmethyl_records(source, "m", None)))

    assert records == [
        ("chr1", 100, 101, "m", ".", 20, 3),
        ("chr1", 200, 201, "m", ".", 20, 18),
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


def test_prepare_bedmethyl_fail_closed_if_output_exists(tmp_path, monkeypatch):
    """Verifies that prepare_bedmethyl fails if output directory already exists."""
    output = tmp_path / "prepared"
    output.mkdir()

    hp1 = tmp_path / "hp1.bedmethyl.gz"
    hp2 = tmp_path / "hp2.bedmethyl.gz"
    ungrouped = tmp_path / "ungrouped.bedmethyl.gz"
    for p in (hp1, hp2, ungrouped):
        write_bedmethyl(p, [bedmethyl_row(100, chrom="chr1")])

    sample_id_b64 = base64.b64encode(b"sample1").decode("utf-8")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(MODULE_PATH),
            "--sample-id-b64", sample_id_b64,
            "--cohort", "control",
            "--task-key", "S000001",
            "--hp1", str(hp1),
            "--hp2", str(hp2),
            "--ungrouped", str(ungrouped),
            "--output", str(output),
        ],
    )

    with pytest.raises(FileExistsError):
        MODULE.main()
