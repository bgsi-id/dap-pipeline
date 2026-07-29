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
vep/
spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz      # optional
spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz.tbi  # optional
```

`reference_release` must change whenever any mounted reference content changes.
It participates in Nextflow task hashes and is written to provenance.

## Required parameters

- `vcf_uri`
- `vcf_index_uri`
- `sample_id`
- `reference_release`
- `output_dir`

Optional SHA-256 parameters are `vcf_sha256` and `vcf_index_sha256`.
`af_cutoff` defaults to `0.01`.

Before normalization, the workflow retains only records on contigs present in
the configured FASTA index. This removes ALT contigs that are incompatible with
the default no-alt GRCh38 reference. The resulting total, retained, and removed
variant counts are published as `<sample_id>.contig-filter.tsv` and included in
provenance.

Keep real input locations in `params.private.json`. This filename is ignored
repository-wide and must never be committed.
