#!/usr/bin/env python3
"""
Validate one sample's phased modkit outputs and create normalized tables.

This module processes per-sample phased bedMethyl files (hp1, hp2, ungrouped),
performs exact-coordinate strand coalescing and deduplication, and creates a
combined all-read dataset using a memory-efficient K-way streaming merge.
"""

import argparse
import base64
import gzip
import hashlib
import heapq
import json
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


PARTITIONS: Tuple[str, ...] = ("hp1", "hp2", "ungrouped")
HEADER: str = "chrom\tstart\tend\tmod_code\tstrand\tvalid_coverage\tmodified\n"

# bedMethyl record tuple: (chrom, start, end, mod_code, strand, valid_coverage, modified_count)
BedMethylRecord = Tuple[str, int, int, str, str, int, int]
# Total ordering key: ((contig_rank, bare_contig), start, end, mod_code, strand)
OrderKey = Tuple[Tuple[int, str], int, int, str, str]
# Coordinate-only key: ((contig_rank, bare_contig), start, end)
CoordinateKey = Tuple[Tuple[int, str], int, int]


def contig_key(name: str) -> Tuple[int, str]:
    """
    Sort contigs numerically by bare name (e.g., 'chr2' -> 2), placing non-standard
    contigs at the end to ensure deterministic coordinate sorting.

    Args:
        name: Chromosome/contig name (e.g., 'chr1', 'chrX', 'chrM').

    Returns:
        A tuple of (rank, bare_name) where standard chromosomes 1-22, X, Y, M/MT
        are ranked 1-25, and non-canonical contigs are assigned rank 1000.
    """
    bare = name.removeprefix("chr")
    order = {str(value): value for value in range(1, 23)}
    order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    return (order.get(bare, 1000), bare)


def order_key(record: BedMethylRecord) -> OrderKey:
    """
    Generate a total ordering key for a bedMethyl record, resolving strand/mod_code ties.

    Args:
        record: A 7-element bedMethyl record tuple.

    Returns:
        A composite key for deterministic sorting.
    """
    chrom, start, end, code, strand, _, _ = record
    return (contig_key(chrom), start, end, code, strand)


def coordinate_key(record: BedMethylRecord) -> CoordinateKey:
    """
    Generate an ordering key based solely on genomic coordinates.

    Args:
        record: A BedMethylRecord tuple containing (chrom, start, end, ...).

    Returns:
        A coordinate-only ordering key ((rank, contig), start, end).
    """
    chrom, start, end = record[:3]
    return (contig_key(chrom), start, end)


def verify_sha256(path: Path, expected: str) -> None:
    """
    Verify a file's SHA-256 digest against an expected manifest value.

    Args:
        path: Path to the target file.
        expected: Expected hex-encoded SHA-256 string (empty string skips check).

    Raises:
        SystemExit: If the calculated digest does not match the expected digest.
    """
    if not expected:
        return
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest().lower() != expected.lower():
        raise SystemExit(f"SHA-256 mismatch for {path.name}")


def parse_region(value: str) -> Optional[Tuple[str, int, int]]:
    """
    Parse a 1-based inclusive genomic interval string (e.g., 'chr22:1-100')
    into 0-based exclusive coordinates (e.g., ('chr22', 0, 100)).

    Args:
        value: Region string in 'contig:start-end' format.

    Returns:
        A tuple of (chrom, 0-based_start, exclusive_end), or None if empty.
    """
    if not value:
        return None
    match = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", value)
    if not match:
        raise SystemExit("region must use contig:start-end with 1-based inclusive coordinates")
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if start < 1 or end < start:
        raise SystemExit("region coordinates are invalid")
    return chrom, start - 1, end


def same_contig(left: str, right: str) -> bool:
    """Check if two contig strings refer to the same chromosome ignoring 'chr' prefix."""
    return left.removeprefix("chr") == right.removeprefix("chr")


def bedmethyl_records(
    path: Path, selected_code: str, region: Optional[Tuple[str, int, int]]
) -> Iterator[BedMethylRecord]:
    """
    Stream a bgzipped 18-column modkit bedMethyl file, yielding validated tuples of:
    (chrom, start, end, mod_code, strand, valid_coverage, modified_count).

    Ensures input is strictly coordinate sorted and filters by an optional genomic region.

    Args:
        path: Path to gzipped bedMethyl input.
        selected_code: Modification code filter (e.g. 'm' for 5mC).
        region: Optional 0-based genomic region tuple (chrom, start, end).

    Yields:
        Validated BedMethylRecord tuples.
    """
    previous = None
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 18:
                raise SystemExit(f"{path.name}:{line_number}: expected 18 bedMethyl columns")
            code = fields[3].split(",", 1)[0]
            if code != selected_code:
                continue
            try:
                record: BedMethylRecord = (
                    fields[0], int(fields[1]), int(fields[2]), code, fields[5],
                    int(fields[9]), int(fields[11]),
                )
            except ValueError as exc:
                raise SystemExit(f"{path.name}:{line_number}: invalid numeric field") from exc
            if record[1] < 0 or record[2] <= record[1]:
                raise SystemExit(f"{path.name}:{line_number}: invalid genomic interval")
            if record[5] < 0 or not 0 <= record[6] <= record[5]:
                raise SystemExit(f"{path.name}:{line_number}: invalid modification counts")
            if region and not (
                same_contig(record[0], region[0])
                and record[2] > region[1]
                and record[1] < region[2]
            ):
                continue
            # bedMethyl is coordinate sorted, but rows sharing an interval are
            # not required to be ordered lexically by modification or strand.
            current = coordinate_key(record)
            if previous is not None and current < previous:
                raise SystemExit(f"{path.name} is not coordinate sorted")
            previous = current
            yield record


def coalesced(
    records: Iterable[BedMethylRecord],
) -> Iterator[BedMethylRecord]:
    """
    Combine (sum) coverage and modified counts for identical coordinates and mod code.

    Upstream pipelines (e.g., nf-core/methylong) produce bedMethyl files where strands
    are combined during pileup. This function performs exact-coordinate strand coalescing
    and deduplication for rows with identical (chrom, start, end, mod_code), preparing
    the data for DSS, which models one count observation per CpG coordinate.

    Args:
        records: An iterable of coordinate-sorted BedMethylRecord tuples.

    Yields:
        Strand-coalesced BedMethylRecord tuples with strand set to '.'.
    """
    current_coordinate = None
    grouped: Dict[Tuple[str, int, int, str], Tuple[int, int]] = {}

    def flush() -> Iterator[BedMethylRecord]:
        for key in sorted(grouped):
            coverage, modified = grouped[key]
            # DSS models one count observation per CpG coordinate and sample.
            yield (*key, ".", coverage, modified)

    for record in records:
        coordinate = record[:3]
        if current_coordinate is not None and coordinate != current_coordinate:
            yield from flush()
            grouped = {}
        current_coordinate = coordinate
        key = record[:4]
        coverage, modified = grouped.get(key, (0, 0))
        grouped[key] = (coverage + record[5], modified + record[6])
    if current_coordinate is not None:
        yield from flush()


def write_records(
    path: Path, records: Iterable[BedMethylRecord]
) -> int:
    """
    Stream normalized records directly to a bgzipped output file.

    Args:
        path: Output file path (.tsv.gz).
        records: Iterable of BedMethylRecord tuples to write.

    Returns:
        Total number of records written.
    """
    count = 0
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER)
        for record in records:
            handle.write("\t".join(map(str, record)) + "\n")
            count += 1
    return count


def normalized_records(path: Path) -> Iterator[BedMethylRecord]:
    """
    Stream tuples from an already-normalized intermediary bedMethyl file.

    Args:
        path: Path to normalized gzipped TSV file.

    Yields:
        BedMethylRecord tuples.
    """
    with gzip.open(path, "rt") as handle:
        next(handle)  # Skip header
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            yield (
                fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4],
                int(fields[5]), int(fields[6]),
            )


def combined_records(
    paths: Sequence[Path],
) -> Iterator[BedMethylRecord]:
    """
    Perform a K-way merge using heapq across multiple partitioned bedMethyl files
    (e.g., hp1, hp2, ungrouped).

    Since each input file is already sorted, this streams across the k inputs with
    O(k) auxiliary memory (k=len(paths)), combining counts for matching coordinates
    to construct a 'primary all-read' combined output.

    Args:
        paths: Sequence of paths to normalized gzipped TSV files.

    Yields:
        Merged BedMethylRecord tuples with aggregated valid_coverage and modified counts.
    """
    iterators = [iter(normalized_records(path)) for path in paths]
    heap: List[Tuple[OrderKey, int, BedMethylRecord]] = []
    for index, iterator in enumerate(iterators):
        try:
            record = next(iterator)
            heapq.heappush(heap, (order_key(record), index, record))
        except StopIteration:
            pass
    while heap:
        key, index, record = heapq.heappop(heap)
        records_at_site = [record]
        try:
            following = next(iterators[index])
            heapq.heappush(heap, (order_key(following), index, following))
        except StopIteration:
            pass
        while heap and heap[0][0] == key:
            _, other_index, other = heapq.heappop(heap)
            records_at_site.append(other)
            try:
                following = next(iterators[other_index])
                heapq.heappush(heap, (order_key(following), other_index, following))
            except StopIteration:
                pass
        yield (
            *record[:5],
            sum(item[5] for item in records_at_site),
            sum(item[6] for item in records_at_site),
        )


def main() -> None:
    """CLI entry point: validate inputs, prepare partitioned files, and create combined output."""
    parser = argparse.ArgumentParser(
        description="Validate one sample's phased modkit outputs and create normalized tables."
    )
    parser.add_argument("--sample-id-b64", required=True, help="Base64-encoded sample ID.")
    parser.add_argument("--cohort", choices=("control", "patient"), required=True,
                        help="Cohort assignment.")
    parser.add_argument("--task-key", required=True, help="Task identification key.")
    parser.add_argument("--hp1", type=Path, required=True, help="Path to HP1 bedMethyl gz.")
    parser.add_argument("--hp2", type=Path, required=True, help="Path to HP2 bedMethyl gz.")
    parser.add_argument("--ungrouped", type=Path, required=True, help="Path to ungrouped bedMethyl gz.")
    parser.add_argument("--hp1-sha256", default="", help="Expected SHA-256 for HP1.")
    parser.add_argument("--hp2-sha256", default="", help="Expected SHA-256 for HP2.")
    parser.add_argument("--ungrouped-sha256", default="", help="Expected SHA-256 for ungrouped.")
    parser.add_argument("--mod-code", default="m", help="Modification code to retain (default: 'm').")
    parser.add_argument("--region", default="", help="Genomic region filter ('contig:start-end').")
    parser.add_argument("--output", type=Path, required=True, help="Output destination directory.")
    args = parser.parse_args()
    region = parse_region(args.region)

    sample_id = base64.b64decode(args.sample_id_b64).decode("utf-8")
    if not sample_id or "\t" in sample_id or "\n" in sample_id:
        raise SystemExit("sample ID is empty or contains a tab/newline")
    # Fail-closed: fail if output directory already exists
    args.output.mkdir()
    source_paths = {name: getattr(args, name) for name in PARTITIONS}
    counts: Dict[str, int] = {}
    normalized_paths: List[Path] = []
    for partition, source in source_paths.items():
        verify_sha256(source, getattr(args, f"{partition}_sha256"))
        output = args.output / f"{partition}.tsv.gz"
        counts[partition] = write_records(
            output, coalesced(bedmethyl_records(source, args.mod_code, region))
        )
        if counts[partition] == 0 and region is None:
            raise SystemExit(f"{source.name} contains no {args.mod_code!r} records")
        normalized_paths.append(output)
    counts["combined"] = write_records(
        args.output / "combined.tsv.gz", combined_records(normalized_paths)
    )
    (args.output / "metadata.json").write_text(
        json.dumps(
            {
                "sample_id": sample_id,
                "cohort": args.cohort,
                "task_key": args.task_key,
                "mod_code": args.mod_code,
                "genomic_region": args.region or "genome-wide",
                "records": counts,
                "sources": {name: path.name for name, path in source_paths.items()},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
