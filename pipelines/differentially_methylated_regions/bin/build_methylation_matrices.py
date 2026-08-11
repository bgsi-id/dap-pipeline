#!/usr/bin/env python3
"""Stream normalized sample tables into coverage-filtered cohort matrices."""

import argparse
import csv
import gzip
import heapq
import json
import math
from pathlib import Path


PARTITIONS = ("combined", "hp1", "hp2", "ungrouped")


def contig_key(name):
    bare = name.removeprefix("chr")
    order = {str(value): value for value in range(1, 23)}
    order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    return (order.get(bare, 1000), bare)


def records(path):
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (
                contig_key(row["chrom"]), int(row["start"]), int(row["end"]),
                row["mod_code"], row["strand"], row["chrom"],
            )
            yield key, int(row["valid_coverage"]), int(row["modified"])


def merged_sites(paths):
    iterators = [iter(records(path)) for path in paths]
    heap = []
    for index, iterator in enumerate(iterators):
        try:
            key, coverage, modified = next(iterator)
            heapq.heappush(heap, (key, index, coverage, modified))
        except StopIteration:
            pass
    while heap:
        key, index, coverage, modified = heapq.heappop(heap)
        values = {index: (coverage, modified)}
        try:
            following = next(iterators[index])
            heapq.heappush(heap, (following[0], index, following[1], following[2]))
        except StopIteration:
            pass
        while heap and heap[0][0] == key:
            _, other_index, other_coverage, other_modified = heapq.heappop(heap)
            values[other_index] = (other_coverage, other_modified)
            try:
                following = next(iterators[other_index])
                heapq.heappush(heap, (following[0], other_index, following[1], following[2]))
            except StopIteration:
                pass
        yield key, values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-coverage", type=int, default=5)
    parser.add_argument("--min-sample-fraction", type=float, default=0.8)
    parser.add_argument("--min-samples-per-cohort", type=int, default=3)
    args = parser.parse_args()
    if args.min_coverage < 1 or not 0 < args.min_sample_fraction <= 1:
        raise SystemExit("invalid coverage or sample-fraction threshold")

    samples = []
    for directory in args.inputs:
        metadata = json.loads((directory / "metadata.json").read_text())
        metadata["directory"] = directory
        samples.append(metadata)
    samples.sort(key=lambda item: (item["cohort"], item["sample_id"]))
    sample_ids = [item["sample_id"] for item in samples]
    if len(sample_ids) != len(set(sample_ids)):
        raise SystemExit("sample IDs must be unique across cohorts")
    cohorts = [item["cohort"] for item in samples]
    regions = {item.get("genomic_region", "genome-wide") for item in samples}
    if len(regions) != 1:
        raise SystemExit("prepared samples disagree on genomic region")
    genomic_region = regions.pop()
    for cohort in ("control", "patient"):
        if cohorts.count(cohort) < args.min_samples_per_cohort:
            raise SystemExit(
                f"{cohort} requires at least {args.min_samples_per_cohort} samples"
            )

    args.output.mkdir()
    metrics = []
    for partition in PARTITIONS:
        destination = args.output / partition
        destination.mkdir()
        with (destination / "metadata.tsv").open("w") as handle:
            handle.write("sample_id\tgroup\n")
            for sample in samples:
                handle.write(f"{sample['sample_id']}\t{sample['cohort']}\n")
        with (destination / "analysis.tsv").open("w") as handle:
            handle.write(f"genomic_region\t{genomic_region}\n")
        paths = [item["directory"] / f"{partition}.tsv.gz" for item in samples]
        retained = observed = 0
        required = {
            cohort: math.ceil(cohorts.count(cohort) * args.min_sample_fraction)
            for cohort in ("control", "patient")
        }
        with gzip.open(destination / "beta.tsv.gz", "wt") as beta, gzip.open(
            destination / "coverage.tsv.gz", "wt"
        ) as coverage_matrix, gzip.open(
            destination / "modified.tsv.gz", "wt"
        ) as modified_matrix, gzip.open(destination / "probes.tsv.gz", "wt") as probes:
            matrix_header = "probe_id\t" + "\t".join(sample_ids) + "\n"
            beta.write(matrix_header)
            coverage_matrix.write(matrix_header)
            modified_matrix.write(matrix_header)
            probes.write("probe_id\tchrom\tposition\n")
            for key, values in merged_sites(paths):
                observed += 1
                passing = {
                    index: value for index, value in values.items()
                    if value[0] >= args.min_coverage
                }
                if any(
                    sum(index in passing for index, value in enumerate(cohorts) if value == cohort)
                    < required[cohort]
                    for cohort in ("control", "patient")
                ):
                    continue
                _, start, _, code, strand, chrom = key
                probe_id = f"{chrom}:{start}:{code}:{strand}"
                beta_row = []
                coverage_row = []
                modified_row = []
                for index in range(len(samples)):
                    if index not in passing:
                        beta_row.append("NA")
                        coverage_row.append("NA")
                        modified_row.append("NA")
                    else:
                        coverage, modified = passing[index]
                        beta_row.append(f"{modified / coverage:.8g}")
                        coverage_row.append(str(coverage))
                        modified_row.append(str(modified))
                beta.write(probe_id + "\t" + "\t".join(beta_row) + "\n")
                coverage_matrix.write(probe_id + "\t" + "\t".join(coverage_row) + "\n")
                modified_matrix.write(probe_id + "\t" + "\t".join(modified_row) + "\n")
                probes.write(f"{probe_id}\t{chrom}\t{start + 1}\n")
                retained += 1
        metrics.append((partition, observed, retained))
    with (args.output.parent / "matrix-metrics.tsv").open("w") as handle:
        handle.write(f"#genomic_region\t{genomic_region}\n")
        handle.write("partition\tsites_observed\tsites_retained\n")
        for row in metrics:
            handle.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
