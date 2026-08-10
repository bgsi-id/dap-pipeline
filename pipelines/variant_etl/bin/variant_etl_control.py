#!/usr/bin/env python3
"""ClickHouse control plane for variant_etl.

Object transfer is deliberately absent. Nextflow stages and publishes files;
ClickHouse reads the published immutable Parquet prefix with its own identity.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlparse

import httpx


def normalize_s3_uri(uri: str) -> str:
    """Collapse duplicate path separators without altering ``s3://``."""
    parsed = urlparse(uri)
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError("published manifest must be an s3:// URI")
    return parsed._replace(path=re.sub(r"/{2,}", "/", parsed.path)).geturl()


def sql(client: httpx.Client, url: str, database: str, statement: str) -> str:
    response = client.post(
        url, params={"database": database}, content=statement.encode("utf-8")
    )
    response.raise_for_status()
    return response.text


def s3_https_glob(manifest_uri: str, region: str, dataset: str) -> str:
    parsed = urlparse(manifest_uri)
    if parsed.scheme != "s3" or not parsed.netloc:
        raise ValueError("published manifest must be an s3:// URI")
    prefix = parsed.path.rsplit("/", 1)[0].lstrip("/")
    return (
        f"https://{parsed.netloc}.s3.{region}.amazonaws.com/"
        f"{prefix}/{dataset}/*/*.parquet"
    )


def load_variants(args: argparse.Namespace) -> None:
    root = Path(args.local_root)
    manifests = list(root.rglob("manifest.json"))
    if len(manifests) != 1:
        raise RuntimeError(f"expected one manifest below {root}, found {len(manifests)}")
    relative = manifests[0].relative_to(root).as_posix()
    manifest_uri = normalize_s3_uri(
        f"{args.published_root.rstrip('/')}/{relative}"
    )

    # Publishing is performed by nf-amazon. ClickHouse may observe the object
    # a few seconds after the process publish completes, so retry only this
    # read/load boundary.
    from variant_ingest.ingest import load_clickhouse

    last_error = None
    for attempt in range(1, 13):
        try:
            load_clickhouse(args.clickhouse_url, args.database, manifest_uri, args.region)
            last_error = None
            break
        except (httpx.HTTPError, RuntimeError) as exc:
            last_error = exc
            if attempt == 12:
                break
            time.sleep(min(5 * attempt, 30))
    if last_error is not None:
        raise last_error

    call_source = s3_https_glob(manifest_uri, args.region, "alt_call")
    safe_manifest = manifest_uri.replace("'", "''")
    safe_sample = args.sample_id.replace("'", "''")
    safe_hash = args.source_sha256.replace("'", "''")
    ddl = f"""
CREATE TABLE IF NOT EXISTS {args.database}.variant_call
(
  release_id LowCardinality(String), batch_id String, sample_id String,
  variant_id String, alt_dosage Int8, genotype String,
  zygosity LowCardinality(String), phased Bool, dp Nullable(Int32),
  qual Nullable(Float64), filter LowCardinality(String), loaded_at DateTime DEFAULT now()
) ENGINE = ReplacingMergeTree(loaded_at)
ORDER BY (release_id, sample_id, variant_id);
CREATE TABLE IF NOT EXISTS {args.database}.variant_ingestion_ledger
(
  manifest_uri String, sample_id String, source_sha256 FixedString(64),
  loaded_at DateTime DEFAULT now()
) ENGINE = ReplacingMergeTree(loaded_at)
ORDER BY (sample_id, source_sha256);
"""
    with httpx.Client(timeout=httpx.Timeout(600, connect=5), trust_env=False) as client:
        # ClickHouse HTTP accepts one DDL statement per request.
        for statement in [part.strip() for part in ddl.split(";") if part.strip()]:
            sql(client, args.clickhouse_url, args.database, statement)
        already = sql(
            client,
            args.clickhouse_url,
            args.database,
            f"SELECT count() FROM {args.database}.variant_ingestion_ledger "
            f"WHERE sample_id='{safe_sample}' AND source_sha256='{safe_hash}'",
        ).strip()
        if int(already or 0) == 0:
            call_count = int(
                sql(
                    client,
                    args.clickhouse_url,
                    args.database,
                    f"SELECT count() FROM s3('{call_source}', 'Parquet')",
                ).strip()
                or 0
            )
            if call_count:
                sql(
                    client,
                    args.clickhouse_url,
                    args.database,
                    f"""INSERT INTO {args.database}.variant_call
(release_id,batch_id,sample_id,variant_id,alt_dosage,genotype,zygosity,phased,dp,qual,filter)
SELECT release_id,batch_id,sample_id,variant_id,alt_dosage,genotype,zygosity,phased,dp,qual,filter
FROM s3('{call_source}', 'Parquet')""",
                )
            sql(
                client,
                args.clickhouse_url,
                args.database,
                f"INSERT INTO {args.database}.variant_ingestion_ledger "
                f"(manifest_uri,sample_id,source_sha256) VALUES "
                f"('{safe_manifest}','{safe_sample}','{safe_hash}')",
            )

    Path(args.output).write_text(
        json.dumps(
            {
                "sample_id": args.sample_id,
                "source_sha256": args.source_sha256,
                "manifest_uri": manifest_uri,
                "status": "loaded",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def info_map(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    if value == ".":
        return result
    for item in value.split(";"):
        key, separator, raw = item.partition("=")
        result[key] = raw if separator else "true"
    return result


def nullable_float(value: str | None) -> float | None:
    if value in (None, "", ".", "-1"):
        return None
    try:
        return float(value.split(",", 1)[0])
    except ValueError:
        return None


def vcf_rows(path: str):
    csq_fields: list[str] = []
    pattern = re.compile(r"Format: ([^\">]+)")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = pattern.search(line)
                if match:
                    csq_fields = match.group(1).strip().split("|")
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(f"malformed VCF row in {path}")
                yield fields[2], info_map(fields[7]), csq_fields


def unique(values):
    return list(dict.fromkeys(value for value in values if value and value != "."))


def present(value: str | None) -> str | None:
    return None if value in (None, "", ".") else value


def detail_values(raw: str | None, fields: list[str]):
    if not raw or raw == "." or not fields:
        return {
            "hgnc": None, "hgnc_id": None, "vep_consequence": [],
            "vep_impact": None, "vep_sift": None, "vep_polyphen": None,
            "vep_transcript_id": None, "vep_hgvsc": None, "vep_hgvsp": None,
            "vep_lof": None, "vep_lof_filter": None, "vep_lof_flags": None,
            "records": [],
        }
    records = [dict(zip(fields, entry.split("|"))) for entry in raw.split(",")]
    consequences = unique(
        consequence
        for record in records
        for consequence in record.get("Consequence", "").split("&")
    )
    preferred = next((record for record in records if record.get("PICK") == "1"), None)
    preferred = preferred or next(
        (record for record in records if present(record.get("MANE_SELECT"))), None
    )
    preferred = preferred or next(
        (record for record in records if record.get("CANONICAL") in {"1", "YES"}), None
    )
    preferred = preferred or records[0]
    return {
        "hgnc": present(preferred.get("SYMBOL")),
        "hgnc_id": present(preferred.get("HGNC_ID")),
        "vep_consequence": consequences,
        "vep_impact": present(preferred.get("IMPACT")),
        "vep_sift": present(preferred.get("SIFT")),
        "vep_polyphen": present(preferred.get("PolyPhen")),
        "vep_transcript_id": present(preferred.get("Feature")),
        "vep_hgvsc": present(preferred.get("HGVSc")),
        "vep_hgvsp": present(preferred.get("HGVSp")),
        "vep_lof": present(preferred.get("LoF")),
        "vep_lof_filter": present(preferred.get("LoF_filter")),
        "vep_lof_flags": present(preferred.get("LoF_flags")),
        "records": records,
    }


def send_json_rows(client, url, database, table, rows, batch_size=5000):
    count = 0
    batch = []
    for row in rows:
        batch.append(json.dumps(row, separators=(",", ":")))
        if len(batch) >= batch_size:
            sql(client, url, database, f"INSERT INTO {table} FORMAT JSONEachRow\n" + "\n".join(batch))
            count += len(batch)
            batch.clear()
    if batch:
        sql(client, url, database, f"INSERT INTO {table} FORMAT JSONEachRow\n" + "\n".join(batch))
        count += len(batch)
    return count


def load_annotations(args: argparse.Namespace) -> None:
    base_table = f"{args.database}.variant_annotation_base"
    detail_table = f"{args.database}.variant_annotation_detail"
    serving_table = f"{args.database}.variant_annotation"
    completion_table = f"{args.database}.{args.completion_table}"
    ddl = [
        f"""CREATE TABLE IF NOT EXISTS {base_table}
(variant_id String, annotation_pack String, gnomad_af Nullable(Float64),
 gnomad_af_max Nullable(Float64), bcsq Array(String), annotated_at DateTime DEFAULT now())
ENGINE = ReplacingMergeTree(annotated_at) ORDER BY (annotation_pack,variant_id)""",
        f"""CREATE TABLE IF NOT EXISTS {detail_table}
(variant_id String, annotation_pack String, hgnc Nullable(String),
 hgnc_id Nullable(String), vep_consequence Array(String),
 vep_impact Nullable(String), vep_sift Nullable(String), vep_polyphen Nullable(String),
 vep_transcript_id Nullable(String), vep_hgvsc Nullable(String), vep_hgvsp Nullable(String),
 vep_lof Nullable(String), vep_lof_filter Nullable(String), vep_lof_flags Nullable(String),
 raw_csq String,
 annotated_at DateTime DEFAULT now())
ENGINE = ReplacingMergeTree(annotated_at) ORDER BY (annotation_pack,variant_id)""",
        f"""CREATE TABLE IF NOT EXISTS {completion_table}
(variant_id String, annotation_pack String, completed_at DateTime DEFAULT now())
ENGINE = ReplacingMergeTree(completed_at) ORDER BY (annotation_pack,variant_id)""",
    ]
    detail_columns = {
        "hgnc_id": "Nullable(String)",
        "vep_impact": "Nullable(String)",
        "vep_sift": "Nullable(String)",
        "vep_transcript_id": "Nullable(String)",
        "vep_hgvsc": "Nullable(String)",
        "vep_hgvsp": "Nullable(String)",
        "vep_lof": "Nullable(String)",
        "vep_lof_filter": "Nullable(String)",
        "vep_lof_flags": "Nullable(String)",
    }
    serving_columns = dict(detail_columns)
    detail_time = datetime.now(timezone.utc).replace(microsecond=0)
    base_time = detail_time - timedelta(seconds=1)
    base_version = base_time.strftime("%Y-%m-%d %H:%M:%S")
    detail_version = detail_time.strftime("%Y-%m-%d %H:%M:%S")

    def base_rows():
        for variant_id, info, _ in vcf_rows(args.base_vcf):
            if not variant_id or variant_id == ".":
                raise ValueError("site VCF lost its canonical variant ID")
            yield {
                "variant_id": variant_id,
                "annotation_pack": args.annotation_pack,
                "gnomad_af": nullable_float(info.get("gnomad_af")),
                "gnomad_af_max": nullable_float(info.get("gnomad_af_max")),
                "bcsq": [] if info.get("BCSQ") in (None, ".") else info["BCSQ"].split(","),
            }

    def serving_base_rows():
        for row in base_rows():
            yield {
                "variant_id": row["variant_id"], "annotation_pack": args.annotation_pack,
                "hgnc": None, "vep_consequence": [], "clinvar_significance": [],
                "clinvar_trait": [], "clinvar_scv": [], "cadd": None,
                "vep_polyphen": None, "phylop": None, "gnomad_af": row["gnomad_af"],
                "annotated_at": base_version,
            }

    def detail_rows():
        for variant_id, info, fields in vcf_rows(args.detail_vcf):
            values = detail_values(info.get("CSQ"), fields)
            yield {
                "variant_id": variant_id, "annotation_pack": args.annotation_pack,
                **{key: value for key, value in values.items() if key != "records"},
                "raw_csq": json.dumps(values["records"], separators=(",", ":")),
                "gnomad_af": nullable_float(info.get("gnomad_af")),
            }

    with httpx.Client(timeout=httpx.Timeout(600, connect=5), trust_env=False) as client:
        for statement in ddl:
            sql(client, args.clickhouse_url, args.database, statement)
        for column, column_type in detail_columns.items():
            sql(client, args.clickhouse_url, args.database,
                f"ALTER TABLE {detail_table} ADD COLUMN IF NOT EXISTS {column} {column_type}")
        for column, column_type in serving_columns.items():
            sql(client, args.clickhouse_url, args.database,
                f"ALTER TABLE {serving_table} ADD COLUMN IF NOT EXISTS {column} {column_type}")
        base_count = send_json_rows(client, args.clickhouse_url, args.database, base_table, base_rows())
        # This projection makes every base site immediately visible to the API
        # and also acts as the annotation-pack anti-join authority.
        send_json_rows(client, args.clickhouse_url, args.database, serving_table, serving_base_rows())
        detail_count = send_json_rows(
            client, args.clickhouse_url, args.database, detail_table,
            ({key: value for key, value in row.items() if key != "gnomad_af"} for row in detail_rows()),
        )
        send_json_rows(
            client,
            args.clickhouse_url,
            args.database,
            serving_table,
            (
                {
                    "variant_id": row["variant_id"], "annotation_pack": args.annotation_pack,
                    "hgnc": row["hgnc"], "hgnc_id": row["hgnc_id"],
                    "vep_consequence": row["vep_consequence"],
                    "vep_impact": row["vep_impact"], "vep_sift": row["vep_sift"],
                    "vep_transcript_id": row["vep_transcript_id"],
                    "vep_hgvsc": row["vep_hgvsc"], "vep_hgvsp": row["vep_hgvsp"],
                    "vep_lof": row["vep_lof"], "vep_lof_filter": row["vep_lof_filter"],
                    "vep_lof_flags": row["vep_lof_flags"],
                    "clinvar_significance": [], "clinvar_trait": [], "clinvar_scv": [],
                    "cadd": None, "vep_polyphen": row["vep_polyphen"], "phylop": None,
                    "gnomad_af": row["gnomad_af"],
                    "annotated_at": detail_version,
                }
                for row in detail_rows()
            ),
        )
        # Completion is committed last. Exporting novel sites anti-joins this
        # table, so a failed partial load cannot permanently hide a site.
        send_json_rows(
            client,
            args.clickhouse_url,
            args.database,
            completion_table,
            (
                {"variant_id": variant_id, "annotation_pack": args.annotation_pack}
                for variant_id, _, _ in vcf_rows(args.base_vcf)
            ),
        )
    Path(args.output).write_text(
        json.dumps(
            {"annotation_pack": args.annotation_pack, "base_rows": base_count,
             "detail_rows": detail_count, "status": "loaded"},
            indent=2, sort_keys=True,
        ) + "\n"
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(dest="command", required=True)
    load = commands.add_parser("load-variants")
    load.add_argument("--local-root", required=True)
    load.add_argument("--published-root", required=True)
    load.add_argument("--clickhouse-url", required=True)
    load.add_argument("--database", required=True)
    load.add_argument("--region", required=True)
    load.add_argument("--sample-id", required=True)
    load.add_argument("--source-sha256", required=True)
    load.add_argument("--output", required=True)
    load.set_defaults(func=load_variants)

    annotation = commands.add_parser("load-annotations")
    annotation.add_argument("--base-vcf", required=True)
    annotation.add_argument("--detail-vcf", required=True)
    annotation.add_argument("--annotation-pack", required=True)
    annotation.add_argument("--assembly", required=True)
    annotation.add_argument("--clickhouse-url", required=True)
    annotation.add_argument("--database", required=True)
    annotation.add_argument("--completion-table", required=True)
    annotation.add_argument("--output", required=True)
    annotation.set_defaults(func=load_annotations)
    return root


if __name__ == "__main__":
    options = parser().parse_args()
    options.func(options)
