nextflow.enable.dsl = 2

/* Development skeleton: reference preparation and projection are pending. */
params.dap_output_uri = null

process DEVELOPMENT_NOTICE {
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true
    output:
    path 'development-notice.txt'
    script:
    """
    printf '%s\\n' 'Genetic ancestry projection (1KGP) is a development skeleton; projection is not implemented.' > development-notice.txt
    """
}

workflow {
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    DEVELOPMENT_NOTICE()
}
