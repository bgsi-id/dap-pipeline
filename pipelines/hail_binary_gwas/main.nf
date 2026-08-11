nextflow.enable.dsl = 2

/* Development skeleton: descriptor and DAP input contract only. */
params.dap_output_uri = null

process DEVELOPMENT_NOTICE {
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true
    output:
    path 'development-notice.txt'
    script:
    """
    printf '%s\\n' 'GWAS (Hail) — binary traits is a development skeleton; association analysis is not implemented.' > development-notice.txt
    """
}

workflow {
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    DEVELOPMENT_NOTICE()
}
