nextflow.enable.dsl = 2

/*
 * TB-Profiler resistance and lineage profiling
 *
 * governed WGS FASTQ (paired- or single-end) -> tb-profiler profile
 *   -> BAM + VCF + resistance/lineage JSON (CSV/TXT when TB-Profiler emits them)
 *
 * All object-store inputs and published outputs are handled by Nextflow.
 * Processes do not invoke an object-store CLI.
 */

params.dap_input_manifest = null
params.dap_output_uri = null
params.batch_id = null


def safe_component(def value, String label) {
    def text = value?.toString()?.trim()
    if (!text || !(text ==~ /[A-Za-z0-9._-]{1,128}/)) {
        error "${label} must contain only A-Z, a-z, 0-9, dot, underscore or hyphen"
    }
    return text
}


def sample_from_asset(def sample) {
    def assets = (sample.assets ?: []).collectEntries { asset ->
        [(asset.role?.toString()?.toLowerCase()): asset]
    }
    if (!sample.id || !assets.fastq_r1) {
        error 'Each resolved sample requires id and a fastq_r1 asset'
    }
    return [
        sample_id: sample.id.toString(),
        reads_r1_uri: assets.fastq_r1.access_uri.toString(),
        reads_r1_sha256: assets.fastq_r1.sha256?.toString(),
        reads_r2_uri: assets.fastq_r2?.access_uri?.toString(),
        reads_r2_sha256: assets.fastq_r2?.sha256?.toString(),
    ]
}


process TBPROFILER_PROFILE {
    tag "${sample_id}"
    container params.tbprofiler_image
    cpus 8
    memory '8 GB'
    time '4h'
    publishDir "${params.dap_output_uri}/results/${sample_id}", mode: 'copy', overwrite: true, failOnError: true

    input:
    tuple val(sample_id), path(reads_r1), val(reads_r1_sha256), path(reads_r2), val(reads_r2_sha256), val(single_end)

    output:
    tuple val(sample_id), path('bam/*.bam'), emit: bam
    tuple val(sample_id), path('vcf/*.vcf.gz'), emit: vcf
    tuple val(sample_id), path('results/*.json'), emit: json
    tuple val(sample_id), path('results/*.csv'), optional: true, emit: csv
    tuple val(sample_id), path('results/*.txt'), optional: true, emit: txt

    script:
    def input_reads = single_end ? "--read1 '${reads_r1}'" : "--read1 '${reads_r1}' --read2 '${reads_r2}'"
    def verify_r1 = reads_r1_sha256 ? "echo '${reads_r1_sha256}  ${reads_r1}' | sha256sum -c -" : 'true'
    def verify_r2 = (!single_end && reads_r2_sha256) ? "echo '${reads_r2_sha256}  ${reads_r2}' | sha256sum -c -" : 'true'
    """
    set -euo pipefail
    ${verify_r1}
    ${verify_r2}
    tb-profiler profile \
      ${input_reads} \
      --prefix '${sample_id}' \
      --threads ${task.cpus}
    """
}


process COLLECT_RUN_SUMMARY {
    container params.tbprofiler_image
    cpus 1
    memory '1 GB'
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true, failOnError: true

    input:
    val sample_ids

    output:
    path 'results/run.json'

    script:
    def ids = sample_ids.collect { "\"${it}\"" }.join(',')
    """
    set -euo pipefail
    mkdir -p results
    cat > results/run.json <<JSON
{
  "pipeline": "tbprofiler-profile",
  "batch_id": "${params.batch_id}",
  "samples": [${ids}]
}
JSON
    """
}


workflow {
    if (!params.dap_input_manifest) error 'dap_input_manifest is required'
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    if (!params.tbprofiler_image) error 'tbprofiler_image is required'
    safe_component(params.batch_id, 'batch_id')

    def resolved = new groovy.json.JsonSlurper().parse(file(params.dap_input_manifest).toFile())
    if (resolved.schema != 'urn:bgsi:dap:resolved-inputs:1' || !(resolved.samples instanceof List)) {
        error 'dap_input_manifest must use resolved-inputs v1 with a samples list'
    }
    def records = resolved.samples.collect { sample_from_asset(it) }
    if (!records) error 'No input samples were resolved'

    def sampleIds = records.collect { it.sample_id }
    if (sampleIds.any { !it } || sampleIds.toSet().size() != sampleIds.size()) {
        error 'Input sample IDs must be present and unique'
    }

    sample_inputs = Channel.fromList(records).map { record ->
        def sampleId = safe_component(record.sample_id, 'sample_id')
        def r1Uri = record.reads_r1_uri
        if (!r1Uri) error "Sample ${sampleId} requires a fastq_r1 URI"
        def r2Uri = record.reads_r2_uri
        def r1File = file(r1Uri, checkIfExists: true)
        def r2File = r2Uri ? file(r2Uri, checkIfExists: true) : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
        tuple(
            sampleId,
            r1File,
            record.reads_r1_sha256?.toString()?.toLowerCase(),
            r2File,
            record.reads_r2_sha256?.toString()?.toLowerCase(),
            !r2Uri,
        )
    }

    profiled = TBPROFILER_PROFILE(sample_inputs)
    COLLECT_RUN_SUMMARY(profiled.bam.map { sample_id, _bam -> sample_id }.collect())
}
