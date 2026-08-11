#!/usr/bin/env python3
"""Validate one sample's phased modkit outputs and create normalized tables."""

import argparse
import base64
import gzip
import hashlib
import heapq
import json
import re
from pathlib import Path


PARTITIONS = ("hp1", "hp2", "ungrouped")
HEADER = "chrom\tstart\tend\tmod_code\tstrand\tvalid_coverage\tmodified\n"


def contig_key(name):
    bare = name.removeprefix("chr")
    order = {str(value): value for value in range(1, 23)}
    order.update({"X": 23, "Y": 24, "M": 25, "MT": 25})
    return (order.get(bare, 1000), bare)


def order_key(record):
    chrom, start, end, code, strand, _, _ = record
    return (contig_key(chrom), start, end, code, strand)


def coordinate_key(record):
    chrom, start, end = record[:3]
    return (contig_key(chrom), start, end)


def verify_sha256(path, expected):
    if not expected:
        return
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest().lower() != expected.lower():
        raise SystemExit(f"SHA-256 mismatch for {path.name}")


def parse_region(value):
    if not value:
        return None
    match = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", value)
    if not match:
        raise SystemExit("region must use contig:start-end with 1-based inclusive coordinates")
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if start < 1 or end < start:
        raise SystemExit("region coordinates are invalid")
    return chrom, start - 1, end


def same_contig(left, right):
    return left.removeprefix("chr") == right.removeprefix("chr")


def bedmethyl_records(path, selected_code, region):
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
                record = (
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


def coalesced(records):
    current_coordinate = None
    grouped = {}

    def flush():
        for key in sorted(grouped):
            coverage, modified = grouped[key]
            yield (*key, coverage, modified)

    for record in records:
        coordinate = record[:3]
        if current_coordinate is not None and coordinate != current_coordinate:
            yield from flush()
            grouped = {}
        current_coordinate = coordinate
        key = record[:5]
        coverage, modified = grouped.get(key, (0, 0))
        grouped[key] = (coverage + record[5], modified + record[6])
    if current_coordinate is not None:
        yield from flush()


def write_records(path, records):
    count = 0
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER)
        for record in records:
            handle.write("\t".join(map(str, record)) + "\n")
            count += 1
    return count


def normalized_records(path):
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            yield (
                fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4],
                int(fields[5]), int(fields[6]),
            )


def combined_records(paths):
    iterators = [iter(normalized_records(path)) for path in paths]
    heap = []
    for index, iterator in enumerate(iterators):
        try:
            record = next(iterator)
            heapq.heappush(heap, (order_key(record), index, record))
        except StopIteration:
            pass
    while heap:
        key, index, record = heapq.heappop(heap)
        records = [record]
        try:
            following = next(iterators[index])
            heapq.heappush(heap, (order_key(following), index, following))
        except StopIteration:
            pass
        while heap and heap[0][0] == key:
            _, other_index, other = heapq.heappop(heap)
            records.append(other)
            try:
                following = next(iterators[other_index])
                heapq.heappush(heap, (order_key(following), other_index, following))
            except StopIteration:
                pass
        yield (*record[:5], sum(item[5] for item in records), sum(item[6] for item in records))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id-b64", required=True)
    parser.add_argument("--cohort", choices=("control", "patient"), required=True)
    parser.add_argument("--task-key", required=True)
    parser.add_argument("--hp1", type=Path, required=True)
    parser.add_argument("--hp2", type=Path, required=True)
    parser.add_argument("--ungrouped", type=Path, required=True)
    parser.add_argument("--hp1-sha256", default="")
    parser.add_argument("--hp2-sha256", default="")
    parser.add_argument("--ungrouped-sha256", default="")
    parser.add_argument("--mod-code", default="m")
    parser.add_argument("--region", default="")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    region = parse_region(args.region)

    sample_id = base64.b64decode(args.sample_id_b64).decode("utf-8")
    if not sample_id or "\t" in sample_id or "\n" in sample_id:
        raise SystemExit("sample ID is empty or contains a tab/newline")
    args.output.mkdir()
    source_paths = {name: getattr(args, name) for name in PARTITIONS}
    counts = {}
    normalized_paths = []
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
