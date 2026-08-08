# Clinical annotation

Native Nextflow implementation of the single-sample rare-disease annotation
flow in `main.sh`.

## Reference mount contract

The execution environment must mount one immutable reference release read-only
at `/reference` (or set `reference_dir` to another absolute mount path):

```text
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
gnomad_v4.1.zip
clinvar.chr.vcf.gz
clinvar.chr.vcf.gz.tbi
vep/homo_sapiens_merged/116_GRCh38/
# or:
homo_sapiens_merged_vep_116_GRCh38/homo_sapiens_merged/116_GRCh38/

spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz       # optional
spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz.tbi   # optional
<compatible SpliceAI indel VCF>.vcf.gz                         # optional
<compatible SpliceAI indel VCF>.vcf.gz.tbi                     # optional
```

`reference_release` must change whenever any mounted reference content changes.
It participates in Nextflow task hashes and is written to provenance.

`vep_cache_subdir=auto` checks the two layouts shown above and then the reference
root. Set it explicitly when using another relative directory. The precheck
requires the merged-cache hierarchy
`homo_sapiens_merged/116_GRCh38`; merely having a non-empty cache directory is
not sufficient.

The ClinVar VCF must use the FASTA's contig convention and must be split and
left-aligned against this exact FASTA when the immutable reference release is
built:

```bash
bcftools norm \
  -m -any \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -c s \
  -Oz \
  -o clinvar.chr.vcf.gz \
  clinvar.input.chr.vcf.gz
bcftools index -f -t clinvar.chr.vcf.gz
```

Do not normalize ClinVar separately for every sample; reference-build-time
normalization preserves sample throughput.

## Required parameters

- `vcf_uri`
- `vcf_index_uri`
- `sample_id`
- `reference_release`
- `output_dir`

Optional SHA-256 parameters are `vcf_sha256` and `vcf_index_sha256`.
`af_cutoff` defaults to `0.01`; `af_max_cutoff` defaults to `0.02`.
`assay_type` defaults to `wgs`; only WGS runs hard-fail when ClinVar matches
zero input records.

SpliceAI is disabled unless `spliceai_indel_name` is explicitly set. When it is
set, both SNV and indel VCFs and both indexes are mandatory. This avoids
starting the VEP release-116 plugin with its required indel input missing.

`vep_buffer_size` defaults to `5000`. VEP still uses the process CPU allocation
for `--fork`; tune CPU count and buffer size together against observed wall
time and peak memory.

## Execution profiles

Select a container profile unless the execution platform supplies its own
executor and container configuration:

```bash
nextflow run main.nf -profile docker -params-file params.private.json
# or
nextflow run main.nf -profile singularity -params-file params.private.json
```

Both bundled profiles mount `reference_dir` read-only at the same absolute
container path. The configuration also writes Nextflow report, timeline, trace,
and DAG files in the launch directory.

Before normalization, the workflow retains only records on contigs present in
the configured FASTA index. This removes ALT contigs that are incompatible with
the default no-alt GRCh38 reference. The resulting total, retained, and removed
variant counts are published as `<sample_id>.contig-filter.tsv` and included in
provenance.

Additional outputs include:

- `<sample_id>.variant-counts.tsv`: keyed gnomAD, ClinVar, and prefilter metrics;
- `<sample_id>.report-metrics.tsv`: VEP CSQ and final report counts;
- `<sample_id>.vep-metadata.tsv`: resolved cache, VEP tuning, and SpliceAI state;
- `<sample_id>.report.tsv`: includes `CLNSIGCONF`, VEP's `PICK` flag, and a
  `REASON_REPORTED` classification (`IMPACT`, `P_LP`, `CONFLICT`, or a
  combination).

The prefilter intentionally rescues ClinVar conflicts before VEP. The final
report labels them as `CONFLICT` rather than treating the substring
`pathogenicity` as a Pathogenic/Likely-pathogenic classification.

Keep real input locations in `params.private.json`. This filename is ignored
repository-wide and must never be committed.
