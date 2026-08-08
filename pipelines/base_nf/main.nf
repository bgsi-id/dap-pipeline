nextflow.enable.dsl = 2

include { READ_BUNDLES; WRITE_RESULTS } from './modules/bundle'

params.dap_input_manifest = null
params.dap_output_uri = null

process BUNDLE_RECEIPT {
    tag 'bundle-contract'
    cpus 1
    memory '1 GB'
    container 'busybox:1.36.1'

    input:
    val bundle

    output:
    path 'results'

    script:
    def summary = [
        schema: 'urn:bgsi:dap:bundle-receipt:1',
        dataset: bundle.dataset,
        cohorts: bundle.cohorts.collectEntries { name, cohort ->
            [(name): [
                pack_id: cohort.pack_id,
                version_id: cohort.version_id,
                receipt_sha256: cohort.receipt_sha256,
                sample_count: cohort.samples.size(),
            ]]
        },
        release_asset_roles: bundle.release_assets.keySet().sort(),
        manifest_sha256: bundle.manifest_sha256,
    ]
    def receipt = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(summary))
    """
    mkdir -p results
    cat > results/bundle-receipt.json <<'JSON'
${receipt}
JSON
    """
}

workflow {
    if (!params.dap_input_manifest) error 'dap_input_manifest is required'
    if (!params.dap_output_uri) error 'dap_output_uri is required'

    bundles = READ_BUNDLES(Channel.value(file(params.dap_input_manifest, checkIfExists: true)))
    receipt = BUNDLE_RECEIPT(bundles.bundle)
    WRITE_RESULTS(receipt)
}
