#!/usr/bin/env python3
"""
Stream normalized sample tables into coverage-filtered cohort matrices.

This module performs a multi-sample K-way streaming merge across prepared bedMethyl
files for each partition ('combined', 'hp1', 'hp2', 'ungrouped'). It applies coverage
and representation filters per cohort and produces beta methylation, valid coverage,
modified count, and probe coordinate matrices formatted for downstream DSS analysis.
"""

import argparse
import csv
import gzip
import heapq
import json
import math
from pathlib import Path
from typing import Dict, Iterator, List, Sequence, Tuple, TypedDict


PARTITIONS: Tuple[str, ...] = ("combined", "hp1", "hp2", "ungrouped")

# Site key structure: ((contig_rank, contig_name), start, end, mod_code, strand, chrom)
SiteKey = Tuple[Tuple[int, str], int, int, str, str, str]


class SampleMetadata(TypedDict):
    """Metadata for a prepared single-sample methylation dataset."""
    sample_id: str
    cohort: str
    genomic_region: str
    directory: Path


def contig_key(name: str) -> Tuple[int, str]:
    """
    Sort contigs numerically by bare name (e.g., 'chr2' -> 2), placing non-standard
    contigs at the end to ensure deterministic coordinate sorting.

    Args:
        name: Chromosome/contig name (e.g., 'chr1', 'chrX', 'chrUn_gl000214v1').

    Returns:
        A tuple of (rank, bare_name) for sorting.
    """
    bare = name.removeprefix("chr")
    order = {str(value): value for value in range(1, 23)}
    order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    return (order.get(bare, 1000), bare)


def records(path: Path) -> Iterator[Tuple[SiteKey, int, int]]:
    """
    Stream normalized bedMethyl records from a gzipped TSV file.

    Yields:
        Tuples of (site_key, valid_coverage, modified_count) where site_key
        encapsulates full genomic coordinate and modification info for sorting.
    """
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key: SiteKey = (
                contig_key(row["chrom"]),
                int(row["start"]),
                int(row["end"]),
                row["mod_code"],
                row["strand"],
                row["chrom"],
            )
            yield key, int(row["valid_coverage"]), int(row["modified"])


def merged_sites(
    paths: Sequence[Path],
) -> Iterator[Tuple[SiteKey, Dict[int, Tuple[int, int]]]]:
    """
    Perform a K-way streaming merge using heapq across multiple sample files.

    Because each sample file is strictly coordinate-sorted, the min-heap maintains
    the active pointer of each of the k file streams. Pop operations yield records
    at the same genomic site across all samples before advancing the iterators,
    requiring O(k) auxiliary memory where k=len(paths).

    Args:
        paths: Sequence of paths to normalized gzipped sample TSVs.

    Yields:
        Tuples of (site_key, sample_values_map) where sample_values_map maps
        the sample index in `paths` to (valid_coverage, modified_count).
    """
    iterators = [iter(records(path)) for path in paths]
    # Heap stores: (site_key, sample_index, valid_coverage, modified_count)
    heap: List[Tuple[SiteKey, int, int, int]] = []
    for index, iterator in enumerate(iterators):
        try:
            key, coverage, modified = next(iterator)
            heapq.heappush(heap, (key, index, coverage, modified))
        except StopIteration:
            pass
    while heap:
        key, index, coverage, modified = heapq.heappop(heap)
        values: Dict[int, Tuple[int, int]] = {index: (coverage, modified)}
        try:
            following_key, following_cov, following_mod = next(iterators[index])
            heapq.heappush(heap, (following_key, index, following_cov, following_mod))
        except StopIteration:
            pass
        # Drain all records sharing the exact same site key
        while heap and heap[0][0] == key:
            _, other_index, other_coverage, other_modified = heapq.heappop(heap)
            values[other_index] = (other_coverage, other_modified)
            try:
                following_key, following_cov, following_mod = next(iterators[other_index])
                heapq.heappush(heap, (following_key, other_index, following_cov, following_mod))
            except StopIteration:
                pass
        yield key, values


def main() -> None:
    """Parse CLI arguments, validate cohort configurations, and stream matrix files."""
    parser = argparse.ArgumentParser(
        description="Stream normalized sample tables into coverage-filtered cohort matrices."
    )
    parser.add_argument("--inputs", nargs="+", type=Path, required=True,
                        help="Prepared sample directories containing metadata.json and partition TSVs.")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output directory for matrix partitions.")
    parser.add_argument("--min-coverage", type=int, default=5,
                        help="Minimum valid read coverage required per CpG observation (default: 5).")
    parser.add_argument("--min-sample-fraction", type=float, default=0.8,
                        help="Minimum fraction of samples per cohort required to meet min-coverage (default: 0.8).")
    parser.add_argument("--min-samples-per-cohort", type=int, default=3,
                        help="Minimum number of samples required in each cohort (default: 3).")
    args = parser.parse_args()

    if args.min_coverage < 1 or not 0 < args.min_sample_fraction <= 1:
        raise SystemExit("invalid coverage or sample-fraction threshold")

    # Load and validate sample metadata across inputs
    samples: List[SampleMetadata] = []
    for directory in args.inputs:
        raw_meta = json.loads((directory / "metadata.json").read_text())
        sample_meta: SampleMetadata = {
            "sample_id": str(raw_meta["sample_id"]),
            "cohort": str(raw_meta["cohort"]),
            "genomic_region": str(raw_meta.get("genomic_region", "genome-wide")),
            "directory": directory,
        }
        samples.append(sample_meta)
    samples.sort(key=lambda item: (item["cohort"], item["sample_id"]))
    sample_ids = [item["sample_id"] for item in samples]
    if len(sample_ids) != len(set(sample_ids)):
        raise SystemExit("sample IDs must be unique across cohorts")
    cohorts = [item["cohort"] for item in samples]
    regions = {item["genomic_region"] for item in samples}
    if len(regions) != 1:
        raise SystemExit("prepared samples disagree on genomic region")
    genomic_region = regions.pop()
    for cohort in ("control", "patient"):
        if cohorts.count(cohort) < args.min_samples_per_cohort:
            raise SystemExit(
                f"{cohort} requires at least {args.min_samples_per_cohort} samples"
            )

    # Fail-closed: fail if output directory already exists
    args.output.mkdir()
    metrics: List[Tuple[str, int, int]] = []

    # Process each partition independently
    for partition in PARTITIONS:
        destination = args.output / partition
        destination.mkdir()
        with (destination / "metadata.tsv").open("w") as handle:
            handle.write("sample_id\tgroup\n")
            for sample in samples:
                handle.write(f"{sample['sample_id']}\t{sample['cohort']}\n")
        with (destination / "analysis.tsv").open("w") as handle:
            handle.write(f"genomic_region\t{genomic_region}\n")

        paths = [sample["directory"] / f"{partition}.tsv.gz" for sample in samples]
        retained = observed = 0
        # Calculate minimum required samples passing coverage threshold per cohort
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
                # Filter to samples meeting the minimum coverage requirement
                passing = {
                    index: value for index, value in values.items()
                    if value[0] >= args.min_coverage
                }
                # Check if both cohorts have sufficient representation
                if any(
                    sum(index in passing for index, value in enumerate(cohorts) if value == cohort)
                    < required[cohort]
                    for cohort in ("control", "patient")
                ):
                    continue

                _, start, _, code, strand, chrom = key
                probe_id = f"{chrom}:{start}:{code}:{strand}"
                beta_row: List[str] = []
                coverage_row: List[str] = []
                modified_row: List[str] = []

                for index in range(len(samples)):
                    if index not in passing:
                        beta_row.append("NA")
                        coverage_row.append("NA")
                        modified_row.append("NA")
                    else:
                        cov, mod = passing[index]
                        beta_row.append(f"{mod / cov:.8g}")
                        coverage_row.append(str(cov))
                        modified_row.append(str(mod))

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
