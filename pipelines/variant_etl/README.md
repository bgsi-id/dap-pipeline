# Variant ETL

`variant_etl` is the incremental ingestion and site-annotation workflow for
DAP Variant. One workflow invocation accepts one sample, an ordinary batch or
a very large batch; sample extraction is channel-driven and bounded by
`INGEST_SAMPLE.maxForks`.

## Runtime contract

This public pipeline intentionally assumes it is launched by `dap-runtime`.
The pipeline does **not** configure a Kubernetes executor, namespace, service
account, IRSA role, work volume, reference volume or cloud credentials.
`dap-runtime` injects those settings in its platform config and mounts:

- the shared Nextflow work volume at the runtime-selected work directory;
- the GRCh38 reference release read-only at `reference_dir` (default
  `/reference`), including the original VEP merged cache;
- an execution identity that can stage governed inputs and publish only to the
  authorised output prefix.

The deployment injects the trusted private `variant_image` and release ID. They
are not user-visible workflow parameters. The other tool images are public and
pinned in `nextflow.config`.

No process invokes `aws s3`, lists a bucket or embeds a bucket name. nf-amazon
stages VCF inputs and publishes Parquet/VCF outputs. ClickHouse reads an exact
published Parquet prefix using its own workload identity.

## Flow

1. Stage each governed VCF, calculate its SHA-256, verify that staged copy and
   extract single-sample VCF facts to Parquet.
2. Publish immutable, source-hash-addressed Parquet and load `variant_dim` and
   `variant_call` in ClickHouse.
3. After the whole ingestion batch completes, export the union of sites that
   is missing the requested annotation pack.
4. Add gnomAD AF with Echtvar.
5. Add inexpensive site-local `BCSQ` with `bcftools csq --local-csq`.
6. Keep missing-frequency and `AF < af_threshold` sites.
7. Run the original offline VEP only on that reduced set, using the mounted
   merged cache and FASTA. The VEP projection retains every transcript CSQ
   record and derives query columns from PICK, then MANE, then canonical
   transcript priority. Allele number, SIFT and PolyPhen are enabled.
8. Load base and detailed annotation projections into ClickHouse and publish
   the audit VCFs and metrics.

Annotation completeness is committed to a separate table only after every
projection has loaded. Novel-site export anti-joins that completion table, so
a failed partial annotation load does not make a site disappear on a fresh
retry.

`--local-csq` is deliberate: a sites-only VCF has no sample haplotypes. A
future sample-level phased consequence pass must operate on the original
sample calls and is a separate fact from this cohort site annotation.

## Inputs

The normal dap-runtime input is `urn:bgsi:dap:resolved-inputs:1`: one or more
governed samples, each with `vcf` and `vcf_index` assets. No catalog checksum is
required for launch. Nextflow stages the assets, calculates their SHA-256 values
and uses the VCF hash in the deterministic ingestion identity.

Example user-visible parameters are in `params.example.json`; dap-runtime adds
the trusted image, resolved release ID and generated `dap_output_uri`.

## Retry semantics

The extractor run ID is `batch + sample + source SHA-256`. ClickHouse keeps
separate manifest and sample/source ledgers. Re-running the same input can
replace the same published keys but does not create a second logical load.
Nextflow `-resume` improves efficiency; it is not the idempotency authority.

The current detailed predicate is intentionally conservative: gnomAD AF is
missing or below the threshold. ClinVar rescue, high-impact BCSQ rescue and
supplementary VEP plugins can be added to that predicate without changing the
pipeline topology.

LoFTEE columns are reserved in the detailed and serving projections. LoFTEE
itself remains disabled until the mounted reference release contains a plugin
compatible with VEP 116 plus its GRCh38 ancestor FASTA, GERP bigWig and SQL
assets; enabling it without that governed bundle would make non-empty runs
non-reproducible.
