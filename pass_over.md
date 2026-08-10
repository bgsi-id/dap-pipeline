# Variant ETL handover

## Runtime launch contract

`pipelines/variant_etl` is a public Nextflow pipeline and is intended to run
only through `dap-runtime`. It does not configure Kubernetes, IRSA, PVCs,
work directories or cloud credentials itself.

`dap-runtime` must supply:

- the Kubernetes executor, task service account/IRSA and shared workspace PVC;
- the read-only GRCh38 reference mount at `/reference`, including FASTA, FAI,
  GFF3, `gnomad_v4.1.zip`, and the original VEP merged cache;
- nf-amazon configuration and an authorised `dap_output_uri`;
- a deployment-owned `variant_image` value for the DAP Variant ingestion
  image.

`variant_image` and `release_id` are required internally but are not
user-facing GUI parameters. For `pipeline_id=variant-etl`, `dap-runtime`
overwrites both with deployment-owned values after resolving the input pack.
The GUI only supplies ordinary workflow parameters such as:

```json
{
  "batch_id": "TEST",
  "annotation_pack": "grch38-v1",
  "af_threshold": 0.01
}
```

The deployment configures the image with `VARIANT_INGEST_IMAGE`. The server
always overwrites an image or release supplied by a caller, so the private
registry address remains trusted and hidden.

The normal governed input is `urn:bgsi:dap:resolved-inputs:1`. Each resolved
sample has `vcf` and `vcf_index` assets. The pipeline computes content hashes
after nf-amazon stages those assets, so the GUI does not require catalogue
checksums. One, batch and megabatch launches use the same `samples` list.

## ClickHouse pass

The pipeline uses ClickHouse as the serving and incremental-state layer:

1. Each VCF becomes source-hash-addressed Parquet. nf-amazon publishes it;
   ClickHouse reads the exact published Parquet prefix using ClickHouse's own
   workload identity.
2. It loads canonical alleles into `variant_dim` and per-sample facts into
   `variant_call`. Manifest and `(sample_id, source_sha256)` ledgers make
   retries logically idempotent.
3. After the entire batch loads, `annotation-sites` exports the union of
   variants not present in `variant_annotation_complete` for the requested
   `annotation_pack`.
4. Echtvar adds gnomAD AF; `bcftools csq --local-csq` adds inexpensive
   site-local BCSQ. Missing/rare sites then receive original offline VEP.
5. The loader stores base fields in `variant_annotation_base`, detailed VEP
   fields in `variant_annotation_detail`, and the API-facing projection in
   `variant_annotation`. It writes `variant_annotation_complete` only last,
   so a failed partial load is eligible for the next retry.

`variant_annotation_complete` is the novelty authority; do not use the API
projection table for that anti-join, otherwise a partial base write can hide a
site before its VEP detail is available.

Current VEP selection is gnomAD-popmax missing or below `af_threshold`
(default `0.01`). ClinVar rescue, high-impact BCSQ rescue and supplementary
VEP plugins are intentionally deferred; they can extend that predicate without
changing the pipeline topology.
