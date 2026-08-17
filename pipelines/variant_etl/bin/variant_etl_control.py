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
CREATE TABLE IF NOT EXISTS {args.database}.variant_dim
(
  variant_id String, assembly LowCardinality(String), contig LowCardinality(String),
  position Int64, xpos Int64, ref String, alt String, variant_type LowCardinality(String)
) ENGINE = ReplacingMergeTree
ORDER BY variant_id;
CREATE TABLE IF NOT EXISTS {args.database}.variant_call
(
  release_id LowCardinality(String), batch_id String, sample_id String,
  variant_id String, alt_dosage Int8, genotype String,
  zygosity LowCardinality(String), phased Bool, phase_set Nullable(String),
  dp Nullable(Int32), gq Nullable(Int32), ad_ref Nullable(Int32), ad_alt Nullable(Int32),
  vaf Nullable(Float32), qual Nullable(Float64), filter LowCardinality(String),
  loaded_at DateTime DEFAULT now()
) ENGINE = ReplacingMergeTree(loaded_at)
ORDER BY (release_id, sample_id, variant_id);
CREATE TABLE IF NOT EXISTS {args.database}.variant_ingestion_ledger
(
  manifest_uri String, sample_id String, source_sha256 FixedString(64),
  loaded_at DateTime DEFAULT now()
) ENGINE = ReplacingMergeTree(loaded_at)
ORDER BY (sample_id, source_sha256);
"""
    alter_call_statements = [
        f"ALTER TABLE {args.database}.variant_dim ADD COLUMN IF NOT EXISTS variant_type LowCardinality(String)",
        f"ALTER TABLE {args.database}.variant_call ADD COLUMN IF NOT EXISTS phase_set Nullable(String)",
        f"ALTER TABLE {args.database}.variant_call ADD COLUMN IF NOT EXISTS gq Nullable(Int32)",
        f"ALTER TABLE {args.database}.variant_call ADD COLUMN IF NOT EXISTS ad_ref Nullable(Int32)",
        f"ALTER TABLE {args.database}.variant_call ADD COLUMN IF NOT EXISTS ad_alt Nullable(Int32)",
        f"ALTER TABLE {args.database}.variant_call ADD COLUMN IF NOT EXISTS vaf Nullable(Float32)",
    ]
    with httpx.Client(timeout=httpx.Timeout(600, connect=5), trust_env=False) as client:
        # ClickHouse HTTP accepts one DDL statement per request.
        for statement in [part.strip() for part in ddl.split(";") if part.strip()]:
            sql(client, args.clickhouse_url, args.database, statement)
        for statement in alter_call_statements:
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
(release_id,batch_id,sample_id,variant_id,alt_dosage,genotype,zygosity,phased,phase_set,dp,gq,ad_ref,ad_alt,vaf,qual,filter)
SELECT s.release_id, s.batch_id, s.sample_id, s.variant_id, s.alt_dosage, s.genotype, s.zygosity, s.phased, s.phase_set, s.dp, s.gq, s.ad_ref, s.ad_alt, s.vaf, s.qual, s.filter
FROM s3('{call_source}', 'Parquet') AS s
LEFT ANTI JOIN
(
  SELECT release_id, sample_id, variant_id
  FROM {args.database}.variant_call
  WHERE sample_id = '{safe_sample}'
) AS existing
ON s.release_id = existing.release_id AND s.sample_id = existing.sample_id AND s.variant_id = existing.variant_id""",
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


def nullable_int(value: str | None) -> int | None:
    if value in (None, "", ".", "-1"):
        return None
    try:
        return int(value.split(",", 1)[0])
    except ValueError:
        return None


def clinvar_stars_from_revstat(revstat: str | None) -> int | None:
    if not revstat or revstat == ".":
        return None
    rev = revstat.lower()
    if "practice_guideline" in rev:
        return 4
    if "reviewed_by_expert_panel" in rev:
        return 3
    if "criteria_provided" in rev and "multiple_submitters" in rev and "no_conflicts" in rev:
        return 2
    if "criteria_provided" in rev and not ("no_assertion" in rev or "no_criteria" in rev):
        return 1
    return 0


def parse_clinvar_sig(clnsig: str | None) -> list[str]:
    if not clnsig or clnsig == ".":
        return []
    parts = re.split(r"[,/|]", clnsig)
    return unique(p.strip() for p in parts if p.strip() and p.strip() != ".")


def parse_clinvar_traits(clndn: str | None) -> list[str]:
    if not clndn or clndn == ".":
        return []
    parts = re.split(r"[,|]", clndn)
    return unique(p.strip().replace("_", " ") for p in parts if p.strip() and p.strip() not in (".", "not_specified", "not_provided"))


def parse_spliceai(value: str | None) -> tuple[float | None, float | None, float | None, float | None, float | None]:
    if not value or value == ".":
        return None, None, None, None, None
    max_ag, max_al, max_dg, max_dl, max_ds = None, None, None, None, None
    for item in value.split(","):
        parts = item.split("|")
        if len(parts) >= 6:
            try:
                ag = float(parts[2]) if parts[2] not in ("", ".") else None
                al = float(parts[3]) if parts[3] not in ("", ".") else None
                dg = float(parts[4]) if parts[4] not in ("", ".") else None
                dl = float(parts[5]) if parts[5] not in ("", ".") else None
                valid = [x for x in (ag, al, dg, dl) if x is not None]
                ds = max(valid) if valid else None
                if ds is not None and (max_ds is None or ds > max_ds):
                    max_ag, max_al, max_dg, max_dl, max_ds = ag, al, dg, dl, ds
            except ValueError:
                continue
    return max_ds, max_ag, max_al, max_dg, max_dl


def parse_base_info(info: dict[str, str]):
    ds_max, ds_ag, ds_al, ds_dg, ds_dl = parse_spliceai(info.get("SpliceAI"))
    scores: dict[str, float] = {}
    if ds_max is not None:
        scores["spliceai_ds_max"] = ds_max
    for key in ("revel", "REVEL", "cadd_phred", "CADD_PHRED", "phylop", "phyloP"):
        val = info.get(key)
        fval = nullable_float(val)
        if fval is not None:
            scores[key.lower()] = fval
    attributes: dict[str, str] = {}
    for key in ("CLNSIGCONF", "CLNVI", "CLNREVSTAT"):
        val = info.get(key)
        if val and val != ".":
            attributes[key.lower()] = val
    return {
        # The prepared gnomAD archive exposes the cohort maximum as its primary
        # AF. Preserve the public gnomad_af field while retaining the source name.
        "gnomad_af": nullable_float(info.get("gnomad_af") or info.get("gnomad_af_max")),
        "gnomad_af_max": nullable_float(info.get("gnomad_af_max")),
        "gnomad_af_popmax": nullable_float(
            info.get("gnomad_af_popmax") or info.get("gnomad_popmax_af") or info.get("popmax_af")
        ),
        "gnomad_nhomalt": nullable_int(
            info.get("gnomad_nhomalt") or info.get("gnomad_nhomalts") or info.get("nhomalt") or info.get("nhomalts")
        ),
        "clinvar_sig": parse_clinvar_sig(info.get("CLNSIG")),
        "clinvar_stars": clinvar_stars_from_revstat(info.get("CLNREVSTAT")),
        "clinvar_traits": parse_clinvar_traits(info.get("CLNDN")),
        "spliceai_ds_max": ds_max,
        "spliceai_ag": ds_ag,
        "spliceai_al": ds_al,
        "spliceai_dg": ds_dg,
        "spliceai_dl": ds_dl,
        "bcsq": [] if info.get("BCSQ") in (None, ".") else info["BCSQ"].split(","),
        "scores": scores,
        "attributes": attributes,
    }


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
            "scores": {}, "attributes": {},
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

    scores: dict[str, float] = {}
    for key in ("REVEL", "REVEL_score", "revel"):
        val = preferred.get(key)
        fval = nullable_float(val)
        if fval is not None:
            scores["revel"] = fval
            break
    for key in ("CADD_PHRED", "cadd_phred"):
        val = preferred.get(key)
        fval = nullable_float(val)
        if fval is not None:
            scores["cadd_phred"] = fval
            break
    for key in ("CADD_RAW", "cadd_raw"):
        val = preferred.get(key)
        fval = nullable_float(val)
        if fval is not None:
            scores["cadd_raw"] = fval
            break
    for key in ("phylop", "phyloP", "phyloP100way_vertebrate"):
        val = preferred.get(key)
        fval = nullable_float(val)
        if fval is not None:
            scores["phylop"] = fval
            break
    sift_raw = preferred.get("SIFT")
    if sift_raw:
        m = re.search(r"\(([0-9.]+)\)", sift_raw)
        if m:
            try:
                scores["sift_score"] = float(m.group(1))
            except ValueError:
                pass
    polyphen_raw = preferred.get("PolyPhen")
    if polyphen_raw:
        m = re.search(r"\(([0-9.]+)\)", polyphen_raw)
        if m:
            try:
                scores["polyphen_score"] = float(m.group(1))
            except ValueError:
                pass

    attributes: dict[str, str] = {}
    for key in ("BIOTYPE", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CANONICAL", "EXON", "INTRON",
                "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons",
                "CLNSIGCONF", "CLNVI", "CLNREVSTAT"):
        val = preferred.get(key)
        if val and val != ".":
            attributes[key.lower()] = val

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
        "scores": scores,
        "attributes": attributes,
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
(
  variant_id String,
  annotation_pack String,
  gnomad_af Nullable(Float64),
  gnomad_af_max Nullable(Float64),
  gnomad_af_popmax Nullable(Float64),
  gnomad_nhomalt Nullable(UInt32),
  clinvar_sig Array(String),
  clinvar_stars Nullable(UInt8),
  clinvar_traits Array(String),
  spliceai_ds_max Nullable(Float32),
  spliceai_ag Nullable(Float32),
  spliceai_al Nullable(Float32),
  spliceai_dg Nullable(Float32),
  spliceai_dl Nullable(Float32),
  bcsq Array(String),
  scores Map(String, Float32),
  attributes Map(String, String),
  annotated_at DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(annotated_at) ORDER BY (annotation_pack, variant_id)""",
        f"""CREATE TABLE IF NOT EXISTS {detail_table}
(
  variant_id String,
  annotation_pack String,
  hgnc Nullable(String),
  hgnc_id Nullable(String),
  vep_consequence Array(String),
  vep_impact Nullable(String),
  vep_sift Nullable(String),
  vep_polyphen Nullable(String),
  vep_transcript_id Nullable(String),
  vep_hgvsc Nullable(String),
  vep_hgvsp Nullable(String),
  vep_lof Nullable(String),
  vep_lof_filter Nullable(String),
  vep_lof_flags Nullable(String),
  scores Map(String, Float32),
  attributes Map(String, String),
  transcript_records_json String,
  annotated_at DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(annotated_at) ORDER BY (annotation_pack, variant_id)""",
        f"""CREATE TABLE IF NOT EXISTS {serving_table}
(
  variant_id String,
  annotation_pack String,
  hgnc Nullable(String),
  hgnc_id Nullable(String),
  vep_consequence Array(String),
  vep_impact Nullable(String),
  vep_sift Nullable(String),
  vep_polyphen Nullable(String),
  vep_transcript_id Nullable(String),
  vep_hgvsc Nullable(String),
  vep_hgvsp Nullable(String),
  vep_lof Nullable(String),
  vep_lof_filter Nullable(String),
  vep_lof_flags Nullable(String),
  gnomad_af Nullable(Float64),
  gnomad_af_max Nullable(Float64),
  gnomad_af_popmax Nullable(Float64),
  gnomad_nhomalt Nullable(UInt32),
  clinvar_sig Array(String),
  clinvar_stars Nullable(UInt8),
  clinvar_traits Array(String),
  spliceai_ds_max Nullable(Float32),
  spliceai_ag Nullable(Float32),
  spliceai_al Nullable(Float32),
  spliceai_dg Nullable(Float32),
  spliceai_dl Nullable(Float32),
  scores Map(String, Float32),
  attributes Map(String, String),
  transcript_records_json String,
  annotated_at DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(annotated_at) ORDER BY (annotation_pack, variant_id)""",
        f"""CREATE TABLE IF NOT EXISTS {completion_table}
(
  variant_id String,
  annotation_pack String,
  completed_at DateTime DEFAULT now()
)
ENGINE = ReplacingMergeTree(completed_at) ORDER BY (annotation_pack, variant_id)""",
    ]
    base_columns = {
        "gnomad_af": "Nullable(Float64)",
        "gnomad_af_max": "Nullable(Float64)",
        "gnomad_af_popmax": "Nullable(Float64)",
        "gnomad_nhomalt": "Nullable(UInt32)",
        "clinvar_sig": "Array(String)",
        "clinvar_stars": "Nullable(UInt8)",
        "clinvar_traits": "Array(String)",
        "spliceai_ds_max": "Nullable(Float32)",
        "spliceai_ag": "Nullable(Float32)",
        "spliceai_al": "Nullable(Float32)",
        "spliceai_dg": "Nullable(Float32)",
        "spliceai_dl": "Nullable(Float32)",
        "scores": "Map(String, Float32)",
        "attributes": "Map(String, String)",
    }
    detail_columns = {
        "hgnc": "Nullable(String)",
        "hgnc_id": "Nullable(String)",
        "vep_consequence": "Array(String)",
        "vep_impact": "Nullable(String)",
        "vep_sift": "Nullable(String)",
        "vep_polyphen": "Nullable(String)",
        "vep_transcript_id": "Nullable(String)",
        "vep_hgvsc": "Nullable(String)",
        "vep_hgvsp": "Nullable(String)",
        "vep_lof": "Nullable(String)",
        "vep_lof_filter": "Nullable(String)",
        "vep_lof_flags": "Nullable(String)",
        "scores": "Map(String, Float32)",
        "attributes": "Map(String, String)",
        "transcript_records_json": "String",
    }
    serving_columns = {
        **detail_columns,
        **base_columns,
    }
    detail_time = datetime.now(timezone.utc).replace(microsecond=0)
    base_time = detail_time - timedelta(seconds=1)
    base_version = base_time.strftime("%Y-%m-%d %H:%M:%S")
    detail_version = detail_time.strftime("%Y-%m-%d %H:%M:%S")

    def base_rows():
        for variant_id, info, _ in vcf_rows(args.base_vcf):
            if not variant_id or variant_id == ".":
                raise ValueError("site VCF lost its canonical variant ID")
            parsed = parse_base_info(info)
            yield {
                "variant_id": variant_id,
                "annotation_pack": args.annotation_pack,
                **parsed,
            }

    def serving_base_rows():
        for row in base_rows():
            yield {
                "variant_id": row["variant_id"],
                "annotation_pack": args.annotation_pack,
                "hgnc": None,
                "hgnc_id": None,
                "vep_consequence": [],
                "vep_impact": None,
                "vep_sift": None,
                "vep_polyphen": None,
                "vep_transcript_id": None,
                "vep_hgvsc": None,
                "vep_hgvsp": None,
                "vep_lof": None,
                "vep_lof_filter": None,
                "vep_lof_flags": None,
                "gnomad_af": row["gnomad_af"],
                "gnomad_af_max": row["gnomad_af_max"],
                "gnomad_af_popmax": row["gnomad_af_popmax"],
                "gnomad_nhomalt": row["gnomad_nhomalt"],
                "clinvar_sig": row["clinvar_sig"],
                "clinvar_stars": row["clinvar_stars"],
                "clinvar_traits": row["clinvar_traits"],
                "spliceai_ds_max": row["spliceai_ds_max"],
                "spliceai_ag": row["spliceai_ag"],
                "spliceai_al": row["spliceai_al"],
                "spliceai_dg": row["spliceai_dg"],
                "spliceai_dl": row["spliceai_dl"],
                "scores": row["scores"],
                "attributes": row["attributes"],
                "transcript_records_json": "[]",
                "annotated_at": base_version,
            }

    def detail_rows():
        for variant_id, info, fields in vcf_rows(args.detail_vcf):
            values = detail_values(info.get("CSQ"), fields)
            base_parsed = parse_base_info(info)
            merged_scores = {**base_parsed["scores"], **values["scores"]}
            merged_attributes = {**base_parsed["attributes"], **values["attributes"]}
            yield {
                "variant_id": variant_id,
                "annotation_pack": args.annotation_pack,
                "hgnc": values["hgnc"],
                "hgnc_id": values["hgnc_id"],
                "vep_consequence": values["vep_consequence"],
                "vep_impact": values["vep_impact"],
                "vep_sift": values["vep_sift"],
                "vep_polyphen": values["vep_polyphen"],
                "vep_transcript_id": values["vep_transcript_id"],
                "vep_hgvsc": values["vep_hgvsc"],
                "vep_hgvsp": values["vep_hgvsp"],
                "vep_lof": values["vep_lof"],
                "vep_lof_filter": values["vep_lof_filter"],
                "vep_lof_flags": values["vep_lof_flags"],
                "gnomad_af": base_parsed["gnomad_af"],
                "gnomad_af_max": base_parsed["gnomad_af_max"],
                "gnomad_af_popmax": base_parsed["gnomad_af_popmax"],
                "gnomad_nhomalt": base_parsed["gnomad_nhomalt"],
                "clinvar_sig": base_parsed["clinvar_sig"],
                "clinvar_stars": base_parsed["clinvar_stars"],
                "clinvar_traits": base_parsed["clinvar_traits"],
                "spliceai_ds_max": base_parsed["spliceai_ds_max"],
                "spliceai_ag": base_parsed["spliceai_ag"],
                "spliceai_al": base_parsed["spliceai_al"],
                "spliceai_dg": base_parsed["spliceai_dg"],
                "spliceai_dl": base_parsed["spliceai_dl"],
                "scores": merged_scores,
                "attributes": merged_attributes,
                "transcript_records_json": json.dumps(values["records"], separators=(",", ":")),
                "annotated_at": detail_version,
            }

    with httpx.Client(timeout=httpx.Timeout(600, connect=5), trust_env=False) as client:
        for statement in ddl:
            sql(client, args.clickhouse_url, args.database, statement)
        for column, column_type in base_columns.items():
            sql(client, args.clickhouse_url, args.database,
                f"ALTER TABLE {base_table} ADD COLUMN IF NOT EXISTS {column} {column_type}")
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
            (
                {
                    "variant_id": row["variant_id"],
                    "annotation_pack": row["annotation_pack"],
                    "hgnc": row["hgnc"],
                    "hgnc_id": row["hgnc_id"],
                    "vep_consequence": row["vep_consequence"],
                    "vep_impact": row["vep_impact"],
                    "vep_sift": row["vep_sift"],
                    "vep_polyphen": row["vep_polyphen"],
                    "vep_transcript_id": row["vep_transcript_id"],
                    "vep_hgvsc": row["vep_hgvsc"],
                    "vep_hgvsp": row["vep_hgvsp"],
                    "vep_lof": row["vep_lof"],
                    "vep_lof_filter": row["vep_lof_filter"],
                    "vep_lof_flags": row["vep_lof_flags"],
                    "scores": row["scores"],
                    "attributes": row["attributes"],
                    "transcript_records_json": row["transcript_records_json"],
                }
                for row in detail_rows()
            ),
        )
        send_json_rows(
            client,
            args.clickhouse_url,
            args.database,
            serving_table,
            detail_rows(),
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
