/*
 * DAP's runtime resolves authorization before Nextflow starts.  This module
 * converts that immutable resolved-inputs v2 document into one small, plain
 * Nextflow value.  Downstream pipelines decide which named cohort(s) and
 * release assets they require.
 */
def read_bundle(def manifest_file) {
    def manifest = new groovy.json.JsonSlurper().parse(manifest_file.toFile())
    if (manifest.schema != 'urn:bgsi:dap:resolved-inputs:2') {
        error "READ_BUNDLES requires urn:bgsi:dap:resolved-inputs:2"
    }
    if (!(manifest.dataset instanceof Map) || !manifest.dataset.id || !manifest.dataset.release_id) {
        error 'READ_BUNDLES requires dataset.id and dataset.release_id'
    }
    if (!(manifest.cohorts instanceof Map)) {
        error 'READ_BUNDLES requires named cohorts'
    }

    def cohorts = manifest.cohorts.collectEntries { name, cohort ->
        if (!(cohort instanceof Map) || !(cohort.samples instanceof List)) {
            error "READ_BUNDLES cohort ${name} must contain a samples list"
        }
        def samples = cohort.samples.collect { it as String }
        if (samples.any { !it.trim() } || samples.size() != samples.toSet().size()) {
            error "READ_BUNDLES cohort ${name} has blank or duplicate sample IDs"
        }
        [(name as String): [
            pack_id: cohort.pack_id?.toString(),
            version_id: cohort.version_id?.toString(),
            receipt_sha256: cohort.receipt_sha256?.toString(),
            samples: samples,
        ]]
    }

    def asset_entries = (manifest.release_assets ?: []).collect { asset ->
        def role = asset.role?.toString()?.trim()?.toLowerCase()
        if (!role || !asset.access_uri) {
            error 'READ_BUNDLES release assets require role and access_uri'
        }
        [role, [
            role: role,
            access_uri: asset.access_uri.toString(),
            name: asset.name?.toString(),
            sha256: asset.sha256?.toString(),
        ]]
    }
    if (asset_entries.collect { item -> item[0] }.toSet().size() != asset_entries.size()) {
        error 'READ_BUNDLES has duplicate release asset roles'
    }
    def release_assets = asset_entries.collectEntries { item -> [(item[0]): item[1]] }

    return [
        schema: manifest.schema.toString(),
        manifest_sha256: manifest.manifest_sha256?.toString(),
        dataset: [
            id: manifest.dataset.id.toString(),
            release_id: manifest.dataset.release_id.toString(),
            release_version: manifest.dataset.release_version?.toString(),
            manifest_sha256: manifest.dataset.manifest_sha256?.toString(),
        ],
        cohorts: cohorts,
        release_assets: release_assets,
    ]
}

workflow READ_BUNDLES {
    take:
    manifest

    main:
    bundle = manifest.map { manifest_file -> read_bundle(manifest_file) }

    emit:
    bundle
}

/*
 * Publish exactly one final result directory to the output URI supplied by
 * dap-web.  An absent or failed destination is an execution failure, never a
 * successful run with inaccessible output.
 */
process WRITE_RESULTS {
    tag 'project-results'
    cpus 1
    memory '1 GB'
    container 'python:3.12-slim'
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true, failOnError: true

    input:
    path results

    output:
    path results

    script:
    """
    test -n '${params.dap_output_uri}'
    test -d results
    """
}
