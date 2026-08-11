import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "bin" / "build_methylation_matrices.py"
SPEC = importlib.util.spec_from_file_location("build_methylation_matrices", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_prepared_sample(root, sample_id, cohort, coverage, modified):
    directory = root / sample_id
    directory.mkdir()
    (directory / "metadata.json").write_text(
        json.dumps(
            {
                "sample_id": sample_id,
                "cohort": cohort,
                "genomic_region": "chr1:101-200",
            }
        )
    )
    for partition in MODULE.PARTITIONS:
        with gzip.open(directory / f"{partition}.tsv.gz", "wt") as handle:
            handle.write(
                "chrom\tstart\tend\tmod_code\tstrand\tvalid_coverage\tmodified\n"
            )
            handle.write(f"chr1\t100\t101\tm\t.\t{coverage}\t{modified}\n")
    return directory


def read_gzip_tsv(path):
    with gzip.open(path, "rt") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def test_emits_raw_dss_counts_and_beta_matrix(tmp_path, monkeypatch):
    inputs = [
        write_prepared_sample(tmp_path, "C1", "control", 20, 2),
        write_prepared_sample(tmp_path, "C2", "control", 22, 3),
        write_prepared_sample(tmp_path, "C3", "control", 18, 1),
        write_prepared_sample(tmp_path, "P1", "patient", 20, 18),
        write_prepared_sample(tmp_path, "P2", "patient", 22, 19),
        write_prepared_sample(tmp_path, "P3", "patient", 18, 16),
    ]
    output = tmp_path / "matrices"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--inputs",
            *map(str, inputs),
            "--output",
            str(output),
            "--min-coverage",
            "5",
            "--min-sample-fraction",
            "1",
            "--min-samples-per-cohort",
            "3",
        ],
    )

    MODULE.main()

    combined = output / "combined"
    assert read_gzip_tsv(combined / "coverage.tsv.gz") == [
        ["probe_id", "C1", "C2", "C3", "P1", "P2", "P3"],
        ["chr1:100:m:.", "20", "22", "18", "20", "22", "18"],
    ]
    assert read_gzip_tsv(combined / "modified.tsv.gz") == [
        ["probe_id", "C1", "C2", "C3", "P1", "P2", "P3"],
        ["chr1:100:m:.", "2", "3", "1", "18", "19", "16"],
    ]
    assert read_gzip_tsv(combined / "beta.tsv.gz")[1] == [
        "chr1:100:m:.",
        "0.1",
        "0.13636364",
        "0.055555556",
        "0.9",
        "0.86363636",
        "0.88888889",
    ]
