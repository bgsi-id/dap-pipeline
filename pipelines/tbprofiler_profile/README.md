# TB-Profiler Profiling

`tbprofiler_profile` wraps [TB-Profiler](https://jodyphelan.gitbook.io/tb-profiler/)
(via the [nf-core `tbprofiler/profile`](https://nf-co.re/modules/tbprofiler_profile)
module) as a DAP pipeline. Given whole-genome sequencing reads for a
*Mycobacterium tuberculosis* isolate, TB-Profiler aligns them to the H37Rv
reference, calls variants, and matches those variants against curated
resistance-mutation and lineage-marker databases bundled inside the tool. The
output is, per sample: a resistance/lineage call (JSON, plus optional CSV/TXT
summaries), the alignment (BAM) and the variant calls (VCF) it was derived
from.

One workflow invocation accepts one sample or a batch of many; each sample is
profiled independently so the batch fans out and back in without any
inter-sample coupling.

## Runtime contract

This public pipeline intentionally assumes it is launched by `dap-runtime`.
The pipeline does **not** configure a Kubernetes executor, namespace, service
account, IRSA role or cloud credentials. `dap-runtime` injects those settings
in its platform config and mounts the shared Nextflow work volume.

Unlike `variant_etl`, this pipeline needs no reference-genome or database
mount: TB-Profiler ships its H37Rv reference and resistance/lineage databases
inside its own container image, so there is no `/reference` dependency to
manage. The `tb-profiler_image` is a public, version-pinned biocontainer
(`nextflow.config`), not a private per-deployment image.

No process invokes an object-store CLI directly or embeds a bucket name;
Nextflow stages the governed FASTQ inputs and publishes outputs.

## Flow

1. Resolve the manifest into one record per sample: a required `fastq_r1`
   asset and an optional `fastq_r2` asset. Presence of `fastq_r2` decides
   paired- vs single-end mode.
2. Stage each sample's FASTQ file(s) and verify their SHA-256 against the
   manifest when the manifest supplies one.
3. Run `tb-profiler profile` per sample: align to H37Rv, call variants, and
   report the resistance/lineage call.
4. Publish each sample's BAM, VCF, JSON (and CSV/TXT when TB-Profiler emits
   them) under `results/<sample_id>/`, plus a batch-level `results/run.json`
   summary listing the pipeline id, `batch_id` and processed sample IDs.

## Inputs

The normal `dap-runtime` input is `urn:bgsi:dap:resolved-inputs:1`: one or
more governed samples. Each sample requires a `fastq_r1` asset; a `fastq_r2`
asset is optional and selects paired-end mode when present. No catalog
checksum is required for launch — Nextflow stages the assets and, when the
manifest carries a SHA-256, this pipeline verifies the staged copy against
it before profiling.

Example user-visible parameters are in `params.example.json`; dap-runtime
supplies the resolved `dap_input_manifest` and generated `dap_output_uri`.

## Outputs

- `results/<sample_id>/bam/*.bam` — alignment against H37Rv.
- `results/<sample_id>/vcf/*.vcf.gz` — variant calls against H37Rv.
- `results/<sample_id>/results/*.json` — resistance/lineage call (always
  produced).
- `results/<sample_id>/results/*.csv` / `*.txt` — optional formatted
  summaries, when TB-Profiler emits them.
- `results/run.json` — batch-level summary (pipeline id, `batch_id`, sample
  IDs processed).
