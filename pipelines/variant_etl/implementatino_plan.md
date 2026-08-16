# Variant ETL Fortification Plan

This plan addresses the clinical hypothesis paradigm shift and ACMG evaluation requirements by updating the `variant_etl` pipeline, the `dap-variant` ingestion engine, and the ClickHouse schemas.

## Open Questions
- Do you want me to also update the ClickHouse schemas defined in `deploy/eks-production.yaml` and `deploy/eks-staging.yaml` in this pass, or will you handle deployment configuration separately?
- Should I update `dap-web/app/clinical/clickhouse.py` and `dap-variant/variant_ingest/api.py` to use the new extensible Map columns, or focus strictly on the ingestion/ETL side for now?

## Proposed Changes

We will create a new branch named `agy` in both the `dap2` and `dap-pipeline` repositories.

---

### 1. `dap-variant` (Rust & Python Extraction)

We will enhance the Rust Parquet extractor to capture crucial sample-level quality and phasing facts needed for technical exclusion and compound-heterozygous evaluation.

#### [MODIFY] `dap-variant/src/lib.rs`
- Update `call_schema` to include `gq` (Int32), `ad_ref` (Int32), `ad_alt` (Int32), `vaf` (Float32), and `phase_set` (Utf8).
- Parse the `GQ`, `AD`, and `PS` FORMAT tags from the VCF fields.
- Calculate `vaf` as `ad_alt / dp`.

#### [MODIFY] `dap-variant/variant_ingest/ingest.py`
- Update the PyArrow `ALT_CALL` schema to match the new Rust schema.
- Update `load_clickhouse` SQL to include `gq`, `ad_ref`, `ad_alt`, `vaf`, and `phase_set` in the `variant_call` table ingestion.

---

### 2. `variant_etl_control.py` (ClickHouse DDL & Load)

We will re-architect the ClickHouse tables to use an extensible Map-based schema for annotations, allowing zero-migration additions of future in-silico predictors.

#### [MODIFY] `dap-pipeline/pipelines/variant_etl/bin/variant_etl_control.py`
- Redesign `variant_annotation` (and base/detail) tables:
  - Add `gnomad_af_popmax`, `gnomad_nhomalt`
  - Add `clinvar_sig`, `clinvar_stars`, `clinvar_traits`
  - Add `spliceai_ds_max`, `spliceai_ag`, `spliceai_al`, `spliceai_dg`, `spliceai_dl`
  - Replace flat nullable scoring columns with `scores Map(String, Float32)` (e.g., `scores['revel']`, `scores['cadd_phred']`).
  - Add `attributes Map(String, String)`.
  - Rename `raw_csq` JSON string to `transcript_records_json`.
- Update the `variant_call` table schema to include the new `gq`, `ad_ref`, `ad_alt`, `vaf`, and `phase_set` columns.
- Update JSON row generation in `load_annotations` to map the incoming VCF fields to these new Map structures.

---

### 3. `main.nf` (Nextflow Pipeline Orchestration)

We will inject the necessary clinical annotations (gnomAD Popmax, ClinVar, SpliceAI, REVEL) and update the rare-variant filter logic to ensure known pathogenic variants are not excluded.

#### [MODIFY] `dap-pipeline/pipelines/variant_etl/main.nf`
- Add reference parameters for the new assets:
  - `clinvar_name = 'clinvar.chr.vcf.gz'`
  - `spliceai_name = 'spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz'`
  - `gnomad_popmax_name = 'gnomad.v3.1.2.echtvar.popmax.v2.zip'`
  - `revel_name = 'revel-v1.3_all_chromosomes.zip'`
- Add `ANNOTATE_CLINVAR` using `bcftools annotate`.
- Add `ANNOTATE_SPLICEAI` using `bcftools annotate`.
- Update `ANNOTATE_AF` to run Echtvar with both `gnomad_v4.1.zip` and the new popmax database.
- Update `SELECT_DETAIL_SITES` to include a ClinVar/SpliceAI rescue clause alongside the frequency threshold (e.g., `< af_threshold OR ClinVar Pathogenic OR SpliceAI > 0.5`).
- Update `ANNOTATE_DETAIL_VEP` to mount and use the REVEL plugin.

## Verification Plan
1. Validate syntax and compilation of the updated Rust code (`cargo check`).
2. Run `pytest` on `dap-web` or `dap-variant` if tests are present to catch obvious integration issues.
3. Push changes to the `agy` branches so they can be reviewed by `codex`.


- `[x]` Create `agy` branch in `dap2` and `dap-pipeline` repositories.
- `[x]` Update `dap-variant/src/lib.rs` schema and parser logic for `GQ`, `AD`, `VAF`, `PS`.
- `[x]` Update `dap-variant/variant_ingest/ingest.py` schema and load SQL.
- `[x]` Update `dap-pipeline/pipelines/variant_etl/bin/variant_etl_control.py` ClickHouse DDL and load logic.
- `[x]` Update `dap-pipeline/pipelines/variant_etl/main.nf` Nextflow processes and parameters.
- `[x]` Verify Rust code compilation and unit test suites.

---

## Codex Review Directive — Required Before Commit or Test Run

Status: **changes requested**. The current working-tree implementation compiles and its
small parser tests pass, but it contains release-blocking integration regressions. Do
not commit or run it against staging ClickHouse until every P0 item below is resolved.
Keep the work scoped to `variant_etl`, `dap-variant`, the required serving consumers,
and explicit schema/deployment migrations; do not include the unrelated dirty DMR or
other untracked files in either commit.

### P0 — Make variant loading single-owner, ordered, and retry-safe

- `variant_ingest.ingest.load_clickhouse()` now inserts `variant_call` rows, while
  `variant_etl_control.load_variants()` already inserts the same Parquet rows. A
  successful run therefore inserts every call twice.
- Worse, `load_clickhouse()` is called before `variant_etl_control` creates or alters
  `variant_call`. A fresh database has no table, and an existing pre-fortification
  table does not yet have the new columns.
- Keep one owner for call loading. Preferred contract:
  `load_clickhouse()` loads only `variant_dim` and its manifest ledger;
  `variant_etl_control.load_variants()` creates/migrates `variant_call`, inserts calls,
  verifies them, then commits the ingestion ledger.
- Preserve retry safety across failures between dimension insert, call insert, and
  ledger commit. A retry must produce exactly one logical call per
  `(release_id, sample_id, variant_id)` and must not rely on background
  `ReplacingMergeTree` merges to hide duplicate physical rows.
- Add a loader test that mocks ClickHouse responses and proves:
  first run inserts calls once, completed retry inserts zero, and failure-before-ledger
  retry cannot leave duplicate calls.

### P0 — Migrate `variant_dim` before using `variant_type`

- The implementation inserts `variant_dim.variant_type`, but the deployed production
  schema in `deploy/eks-production.yaml` has no such column and no runtime
  `ALTER TABLE ... ADD COLUMN` is issued. This fails before annotations start.
- Either add an idempotent migration before the insert and update both staging and
  production bootstrap schemas, or omit `variant_type` from the ClickHouse projection.
- Verify loading against both an old schema and a fresh schema.

### P0 — Prepare REVEL correctly; never pass the vendor ZIP to VEP

- `revel-v1.3_all_chromosomes.zip` is source material, not a VEP plugin database.
  `--plugin REVEL,<zip>` is invalid.
- Prepare the GRCh38 REVEL file once: unzip it, convert CSV to tab-separated data,
  retain valid GRCh38 coordinates, coordinate-sort it, bgzip it, and create its tabix
  index. Store the prepared `.tsv.gz` and `.tbi` in the mounted reference pack.
- Invoke the release-116 plugin with the documented form
  `--plugin REVEL,file=/reference/.../revel_grch38.tsv.gz` and ensure `REVEL.pm` is
  present in the configured plugin directory. See the official Ensembl plugin
  documentation: https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
- Make the prepared file and index mandatory in `PRECHECK_REFERENCES`; do not silently
  skip a clinical annotation promised by the annotation-pack version.
- Add a tiny VEP smoke test containing one known REVEL hit and assert that `REVEL`
  appears in the emitted CSQ header and value.

### P0 — Introduce a new annotation lineage and backfill existing sites

- The pipeline still defaults to `annotation_pack = grch38-v2`. Existing rows in
  `variant_annotation_complete` will therefore anti-join out previously annotated
  variants, so this implementation will annotate only newly observed sites and leave
  the existing cohort on the old schema/content.
- Assign a new immutable annotation-pack ID (for example `grch38-v3`, with an explicit
  reference manifest/version), update deployment defaults consistently, and run a
  controlled backfill of all sites for the new pack. Do not delete or mutate v2.
- Prove that an existing variant receives v3 base/detail annotations and that v2
  remains queryable for reproducibility.

### P0 — Update serving consumers atomically with the schema

- `dap-variant/variant_ingest/api.py` and
  `dap-web/app/clinical/clickhouse.py` still consume legacy fields such as
  `clinvar_significance`, `clinvar_trait`, `clinvar_scv`, `cadd`, and `phylop`.
  The new loader writes `clinvar_sig`, `clinvar_traits`, `scores`, and `attributes`
  instead, so the Interpretation Console would silently lose ClinVar and predictor
  values even when annotation succeeds.
- Choose and document one compatibility approach:
  update all readers to project values from the new columns/maps, or maintain a
  compatibility projection/materialized view that exposes the old API contract.
- Preserve the public response model during the migration, including CADD, REVEL,
  PolyPhen, phyloP, ClinVar significance/traits/review stars and SCV identifiers.
  `CLNVI` is not an SCV substitute; capture actual assertion/accession identifiers if
  SCV is part of the required clinical contract.
- Add end-to-end fixture coverage from annotated VCF -> ClickHouse JSON row ->
  dap-variant response -> clinical `Finding` mapping.

### P1 — Validate exact reference field aliases

- Do not assume aliases from filenames. Inspect and record the INFO aliases embedded
  in both Echtvar archives. The official Echtvar examples use names such as
  `gnomad_popmax_af` and `gnomad_nhomalts`, while the current filter asks for
  `INFO/gnomad_af_popmax` and the parser does not accept plural `gnomad_nhomalts`.
  An undefined bcftools field can abort selection; a mismatched parser silently makes
  the new columns null.
- Keep multiple `-e` archives: this is supported by Echtvar's official interface
  (https://github.com/brentp/echtvar#annotate). Add a header-contract check immediately
  after annotation and fail with a clear message unless every required alias exists.
- Use one canonical internal name after ingest, independent of the source alias.

### P1 — Make ClinVar and SpliceAI matching/filtering allele-correct

- Require and precheck `.tbi`/`.csi` indexes for ClinVar and SpliceAI, not just their
  data files.
- Confirm that query and annotation VCFs use identical GRCh38 contig naming,
  decomposition, and left-normalization. Add a small known-hit/miss test and report
  match counts; a successful command with zero matches is not success.
- Replace the regex scan of the compound `SpliceAI` string with deterministic numeric
  extraction of `DS_AG`, `DS_AL`, `DS_DG`, and `DS_DL`, then filter on a derived
  `DS_MAX >= threshold` value.
- Define ClinVar rescue semantics explicitly. Preserve benign and conflicting records
  in base annotation, but only rescue detail annotation for the approved set (for
  example P/LP and conflicts containing a P/LP assertion). Do not accidentally rescue
  every conflict merely because `CLNSIGCONF` is present or because the word
  `pathogenicity` appears inside a generic conflict label.

### P1 — Strengthen extraction and schema tests

- The new dap-variant tests currently inspect only the Python fallback parser and
  schema names. Add an end-to-end native Rust extraction test that reads a phased,
  multiallelic VCF and asserts Parquet values for `GQ`, per-ALT `AD`, `VAF`, and `PS`.
- Test missing and malformed FORMAT values and require identical Python/Rust behavior.
- Test the complete `load_annotations` JSON rows and ClickHouse DDL, not only helper
  parsers. Include Map serialization, base-to-detail replacement, empty-detail input,
  and migration from the current deployed schema.

### Required verification evidence

Before marking this plan complete, attach results for all of the following:

1. `cargo check` and native Parquet extraction tests.
2. Full dap-variant tests, including loader idempotency.
3. Variant ETL control tests covering generated DDL/JSON rows and old-schema migration.
4. Nextflow parse/config validation plus a small local or Kubernetes known-variant run.
5. ClickHouse assertions showing no duplicate calls, non-null expected annotation
   fields, and correct v2/v3 lineage counts.
6. API/UI contract tests proving the Interpretation Console receives ClinVar, CADD,
   REVEL, PolyPhen, phyloP, gnomAD popmax, SpliceAI, genotype quality, allele depth,
   VAF, and phase-set values.
7. Commit the intended changes on both `agy` branches. At review time, root `agy`
   still pointed to the same commit as `staging`, and dap-pipeline `agy` still pointed
   to the same commit as `dev`; all reviewed implementation and tests were uncommitted.
