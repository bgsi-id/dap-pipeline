nextflow.enable.dsl = 2

if (!params.vcf_uri || !params.vcf_index_uri) {
    error "vcf_uri and vcf_index_uri are required"
}
if (!params.vcf_sha256 || !params.vcf_index_sha256) {
    error "vcf_sha256 and vcf_index_sha256 are required"
}

process VERIFY_SAMPLE_VCF {
    tag "verify-sample-vcf"

    input:
    path vcf
    path vcf_index

    output:
    path "summary.txt", emit: summary

    script:
    """
    echo '${params.vcf_sha256}  ${vcf}' | sha256sum -c -
    echo '${params.vcf_index_sha256}  ${vcf_index}' | sha256sum -c -
    variants=\$(gzip -cd '${vcf}' | awk '!/^#/ { count += 1 } END { print count + 0 }')
    {
      echo "VCF input verified"
      echo "file=${vcf}"
      echo "index=${vcf_index}"
      echo "variants=\${variants}"
    } | tee summary.txt
    """
}

workflow {
    vcf = Channel.fromPath(params.vcf_uri, checkIfExists: true)
    vcf_index = Channel.fromPath(params.vcf_index_uri, checkIfExists: true)
    VERIFY_SAMPLE_VCF(vcf, vcf_index)
}
