nextflow.enable.dsl = 2

process VERIFY_SAMPLE_VCF {
    tag "verify-${sample_id}"

    input:
    tuple val(sample_id), path(vcf), val(vcf_sha256), path(vcf_index), val(vcf_index_sha256)

    output:
    tuple val(sample_id), path("${sample_id}.summary.txt"), emit: summary

    script:
    """
    echo '${vcf_sha256}  ${vcf}' | sha256sum -c -
    echo '${vcf_index_sha256}  ${vcf_index}' | sha256sum -c -
    variants=\$(gzip -cd '${vcf}' | awk '!/^#/ { count += 1 } END { print count + 0 }')
    {
      echo "VCF input verified"
      echo "file=${vcf}"
      echo "index=${vcf_index}"
      echo "variants=\${variants}"
    } | tee '${sample_id}.summary.txt'
    """
}

workflow {
    if (!params.dap_input_manifest) {
        error "dap_input_manifest is required"
    }

    def manifest = new groovy.json.JsonSlurper().parse(new File(params.dap_input_manifest as String))
    if (manifest.schema != 'urn:bgsi:dap:resolved-inputs:1' || !manifest.samples) {
        error "dap_input_manifest is invalid or empty"
    }

    def resolvedSamples = manifest.samples.collect { sample ->
        def assets = sample.assets.collectEntries { [(it.role): it] }
        if (!assets.vcf || !assets.vcf_index) {
            error "Sample ${sample.id} does not contain vcf and vcf_index assets"
        }
        tuple(
            sample.id as String,
            assets.vcf.access_uri as String,
            assets.vcf.sha256 as String,
            assets.vcf_index.access_uri as String,
            assets.vcf_index.sha256 as String,
        )
    }

    inputs = Channel.fromList(resolvedSamples).map { sample_id, vcf_uri, vcf_sha256, vcf_index_uri, vcf_index_sha256 ->
        tuple(
            sample_id,
            file(vcf_uri, checkIfExists: true),
            vcf_sha256,
            file(vcf_index_uri, checkIfExists: true),
            vcf_index_sha256,
        )
    }
    VERIFY_SAMPLE_VCF(inputs)
}
