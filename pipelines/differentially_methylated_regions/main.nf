nextflow.enable.dsl = 2

include { WRITE_RESULTS } from '../base_nf/modules/bundle'

params.dap_input_manifest = null
params.dap_output_uri = null
params.cohort_id = null
params.mod_code = 'm'
params.genomic_region = ''
params.min_coverage = 5
params.min_sample_fraction = 0.8
params.min_samples_per_cohort = 3
params.dss_p_threshold = 0.001
params.dss_delta = 0.05
params.dss_min_length = 50
params.dss_min_cpgs = 3
params.dss_merge_distance = 50
params.dss_significant_fraction = 0.5
params.dss_equal_dispersion = true
params.dss_smoothing = false
params.dmr_memory = '32 GB'


def safe_analysis_id(def value) {
    def text = value?.toString()?.trim()
    if (!text || !(text ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) {
        error 'cohort_id contains unsupported characters'
    }
    return text
}


def sample_methylation_record(def sample, int index) {
    if (!sample.id || !(sample.cohort in ['control', 'patient'])) {
        error 'Each resolved sample requires id and control/patient cohort membership'
    }
    def candidates = (sample.assets ?: []).findAll {
        it.role?.toString()?.toLowerCase() == 'bedmethyl'
    }
    def suffixes = [
        hp1: '.wf_mods.1.bedmethyl.gz',
        hp2: '.wf_mods.2.bedmethyl.gz',
        ungrouped: '.wf_mods.ungrouped.bedmethyl.gz',
    ]
    def selected = suffixes.collectEntries { partition, suffix ->
        def matches = candidates.findAll { it.name?.toString()?.endsWith(suffix) }
        if (matches.size() != 1) {
            def received = (sample.assets ?: []).collect {
                "${it.role ?: '<missing-role>'}:${it.name ?: '<missing-name>'}"
            }.join(', ')
            error "Sample ${sample.id} requires exactly one bedmethyl asset ending ${suffix}; received [${received}]"
        }
        [(partition): matches.first()]
    }
    selected.each { partition, asset ->
        if (!asset.access_uri) error "Sample ${sample.id} ${partition} asset has no access URI"
        if (asset.sha256 && !(asset.sha256.toString() ==~ /[0-9a-fA-F]{64}/)) {
            error "Sample ${sample.id} ${partition} asset has an invalid SHA-256"
        }
    }
    return [
        task_key: String.format('S%06d', index + 1),
        sample_id: sample.id.toString(),
        cohort: sample.cohort.toString(),
        hp1: selected.hp1,
        hp2: selected.hp2,
        ungrouped: selected.ungrouped,
    ]
}


process PREPARE_SAMPLE_METHYLATION {
    tag "${sample_id}"
    cpus 2
    memory '8 GB'
    container params.python_image

    input:
    tuple val(task_key), val(sample_id), val(cohort), path(hp1), val(hp1_sha256), path(hp2), val(hp2_sha256), path(ungrouped), val(ungrouped_sha256)

    output:
    path "${task_key}.prepared"

    script:
    def encodedSample = sample_id.bytes.encodeBase64().toString()
    """
    set -euo pipefail
    python3 '${projectDir}/bin/prepare_bedmethyl.py' \
      --sample-id-b64 '${encodedSample}' \
      --cohort '${cohort}' \
      --task-key '${task_key}' \
      --hp1 '${hp1}' --hp1-sha256 '${hp1_sha256}' \
      --hp2 '${hp2}' --hp2-sha256 '${hp2_sha256}' \
      --ungrouped '${ungrouped}' --ungrouped-sha256 '${ungrouped_sha256}' \
      --mod-code '${params.mod_code}' \
      --region '${params.genomic_region}' \
      --output '${task_key}.prepared'
    """

    stub:
    """
    mkdir '${task_key}.prepared'
    printf '{"sample_id":"%s","cohort":"%s","task_key":"%s","mod_code":"m"}\n' '${sample_id}' '${cohort}' '${task_key}' > '${task_key}.prepared/metadata.json'
    for partition in combined hp1 hp2 ungrouped; do
      printf 'chrom\tstart\tend\tmod_code\tstrand\tvalid_coverage\tmodified\n' | gzip > "${task_key}.prepared/\${partition}.tsv.gz"
    done
    """
}


process BUILD_METHYLATION_MATRICES {
    tag params.cohort_id
    cpus 2
    memory params.dmr_memory
    container params.python_image

    input:
    path prepared_samples

    output:
    tuple val('combined'), path('matrices/combined'), emit: combined
    tuple val('hp1'), path('matrices/hp1'), emit: hp1
    tuple val('hp2'), path('matrices/hp2'), emit: hp2
    tuple val('ungrouped'), path('matrices/ungrouped'), emit: ungrouped
    path 'matrix-metrics.tsv', emit: metrics

    script:
    def inputs = prepared_samples.collect { "'${it}'" }.join(' ')
    """
    set -euo pipefail
    python3 '${projectDir}/bin/build_methylation_matrices.py' \
      --inputs ${inputs} \
      --output matrices \
      --min-coverage ${params.min_coverage} \
      --min-sample-fraction ${params.min_sample_fraction} \
      --min-samples-per-cohort ${params.min_samples_per_cohort}
    """

    stub:
    """
    mkdir -p matrices/{combined,hp1,hp2,ungrouped}
    for partition in combined hp1 hp2 ungrouped; do
      printf 'probe_id\tS1\n' | gzip > "matrices/\${partition}/beta.tsv.gz"
      printf 'probe_id\tS1\n' | gzip > "matrices/\${partition}/coverage.tsv.gz"
      printf 'probe_id\tS1\n' | gzip > "matrices/\${partition}/modified.tsv.gz"
      printf 'probe_id\tchrom\tposition\n' | gzip > "matrices/\${partition}/probes.tsv.gz"
      printf 'sample_id\tgroup\nS1\tcontrol\n' > "matrices/\${partition}/metadata.tsv"
      printf 'genomic_region\t%s\n' '${params.genomic_region ?: 'genome-wide'}' > "matrices/\${partition}/analysis.tsv"
    done
    printf 'partition\tsites_observed\tsites_retained\n' > matrix-metrics.tsv
    """
}


process CALL_DMRS {
    tag "${params.cohort_id}:${partition}"
    cpus 4
    memory params.dmr_memory
    container params.r_image

    input:
    tuple val(partition), path(matrix_dir)

    output:
    tuple val(partition), path("${partition}.dmr")

    script:
    """
    set -euo pipefail
    Rscript '${projectDir}/bin/call_dmrs.R' \
      '${matrix_dir}' '${partition}.dmr' '${partition}' '${params.genomic_region}' \
      ${params.dss_p_threshold} ${params.dss_delta} ${params.dss_min_length} \
      ${params.dss_min_cpgs} ${params.dss_merge_distance} ${params.dss_significant_fraction} \
      ${params.dss_equal_dispersion} ${params.dss_smoothing} ${task.cpus}
    """

    stub:
    """
    mkdir '${partition}.dmr'
    printf 'chr\tstart\tend\tlength\tnCG\tmeanMethy1\tmeanMethy2\tdiff.Methy\tareaStat\n' > '${partition}.dmr/dmr-results.tsv'
    printf 'chr\tpos\tmu1\tmu2\tdiff\tdiff.se\tstat\tphi1\tphi2\tpval\tfdr\n' > '${partition}.dmr/cpg-results.tsv'
    printf 'chr\tpos\tmu1\tmu2\tdiff\tdiff.se\tstat\tphi1\tphi2\tpval\tfdr\n' > '${partition}.dmr/dml-results.tsv'
    printf 'metric\tvalue\npartition\t${partition}\n' > '${partition}.dmr/dmr-metrics.tsv'
    printf 'partition\tchrom\tposition\tgroup\tmean_beta\n' > '${partition}.dmr/top-dmr-profile.tsv'
    printf 'partition\t${partition}\n' > '${partition}.dmr/provenance.tsv'
    touch '${partition}.dmr/top-dmr-profile.png'
    touch '${partition}.dmr/top-dmr-dss.png'
    """
}


process COLLECT_RESULTS {
    tag params.cohort_id
    cpus 1
    memory '2 GB'
    container params.python_image

    input:
    path partition_results
    path matrix_metrics

    output:
    path 'results'

    script:
    """
    set -euo pipefail
    mkdir results
    cp '${matrix_metrics}' results/${params.cohort_id}.matrix-metrics.tsv
    for source in ${partition_results.collect { "'${it}'" }.join(' ')}; do
      partition="\${source%.dmr}"
      mkdir "results/\${partition}"
      cp "\${source}"/* "results/\${partition}/"
    done
    """
}


workflow {
    if (!params.dap_input_manifest) error 'dap_input_manifest is required'
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    safe_analysis_id(params.cohort_id)
    if (!(params.mod_code ==~ /[A-Za-z][A-Za-z0-9?+-]{0,15}/)) error 'invalid mod_code'
    if (params.genomic_region && !(params.genomic_region ==~ /[A-Za-z0-9_.-]+:[0-9]+-[0-9]+/)) error 'genomic_region must use contig:start-end'
    if ((params.min_coverage as Integer) < 1) error 'min_coverage must be positive'
    if ((params.min_samples_per_cohort as Integer) < 3) error 'min_samples_per_cohort must be at least 3'
    if ((params.min_sample_fraction as BigDecimal) <= 0 || (params.min_sample_fraction as BigDecimal) > 1) error 'min_sample_fraction must be in (0,1]'
    if ((params.dss_p_threshold as BigDecimal) <= 0 || (params.dss_p_threshold as BigDecimal) >= 1) error 'dss_p_threshold must be in (0,1)'
    if ((params.dss_delta as BigDecimal) < 0 || (params.dss_delta as BigDecimal) >= 1) error 'dss_delta must be in [0,1)'
    if ((params.dss_min_length as Integer) < 1) error 'dss_min_length must be positive'
    if ((params.dss_min_cpgs as Integer) < 2) error 'dss_min_cpgs must be at least 2'
    if ((params.dss_merge_distance as Integer) < 0) error 'dss_merge_distance must be non-negative'
    if ((params.dss_significant_fraction as BigDecimal) <= 0 || (params.dss_significant_fraction as BigDecimal) > 1) error 'dss_significant_fraction must be in (0,1]'

    def resolved = new groovy.json.JsonSlurper().parse(file(params.dap_input_manifest).toFile())
    if (resolved.schema != 'urn:bgsi:dap:resolved-inputs:2' || !(resolved.samples instanceof List)) {
        error 'DMR requires a resolved-inputs v2 manifest with governed sample assets'
    }
    if (!(resolved.cohorts?.control?.samples instanceof List) || !(resolved.cohorts?.patient?.samples instanceof List)) {
        error 'DMR requires control and patient cohort packs'
    }
    def records = resolved.samples.withIndex().collect { sample, index ->
        sample_methylation_record(sample, index)
    }
    if (!records) error 'No methylation samples were resolved'
    def ids = records.collect { it.sample_id }
    if (ids.toSet().size() != ids.size()) error 'Resolved sample IDs must be unique'
    ['control', 'patient'].each { cohort ->
        def expected = (resolved.cohorts[cohort].samples as List).collect { it.toString() }.toSet()
        def observed = records.findAll { it.cohort == cohort }.collect { it.sample_id }.toSet()
        if (expected != observed) error "Resolved ${cohort} membership and sample assets disagree"
    }

    sample_inputs = Channel.fromList(records).map { record ->
        tuple(
            record.task_key,
            record.sample_id,
            record.cohort,
            file(record.hp1.access_uri.toString(), checkIfExists: true),
            record.hp1.sha256?.toString() ?: '',
            file(record.hp2.access_uri.toString(), checkIfExists: true),
            record.hp2.sha256?.toString() ?: '',
            file(record.ungrouped.access_uri.toString(), checkIfExists: true),
            record.ungrouped.sha256?.toString() ?: '',
        )
    }
    prepared = PREPARE_SAMPLE_METHYLATION(sample_inputs)
    matrices = BUILD_METHYLATION_MATRICES(prepared.collect())
    partition_matrices = matrices.combined.mix(matrices.hp1, matrices.hp2, matrices.ungrouped)
    partition_results = CALL_DMRS(partition_matrices)
    results = COLLECT_RESULTS(partition_results.map { partition, directory -> directory }.collect(), matrices.metrics)
    WRITE_RESULTS(results)
}
