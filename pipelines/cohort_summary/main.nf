nextflow.enable.dsl = 2

params.dap_input_manifest = null
params.dap_output_uri = null

process WRITE_COHORT_SUMMARY {
    container params.python_image
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true

    input:
    path manifest

    output:
    path 'cohort-summary.tsv'
    path 'cohort-summary.json'

    script:
    """
    python3 '${projectDir}/bin/summarize_cohort.py' '${manifest}' cohort-summary.tsv cohort-summary.json
    """
}

workflow {
    if (!params.dap_input_manifest) error 'dap_input_manifest is required'
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    WRITE_COHORT_SUMMARY(file(params.dap_input_manifest, checkIfExists: true))
}
