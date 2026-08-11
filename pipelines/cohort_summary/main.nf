nextflow.enable.dsl = 2

params.dap_input_manifest = null
params.dap_output_uri = null

process WRITE_COHORT_SUMMARY {
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true

    input:
    path manifest

    output:
    path 'cohort-summary.tsv'

    script:
    """
    python3 -c 'import json; data=json.load(open("${manifest}")); samples=data.get("samples") or []; rows=["sample_id\\tasset_roles"]; rows += ["{}\\t{}".format(sample.get("id", ""), ",".join(sorted(asset.get("role", "") for asset in sample.get("assets", [])))) for sample in samples]; open("cohort-summary.tsv", "w").write("\\n".join(rows) + "\\n")'
    """
}

workflow {
    if (!params.dap_input_manifest) error 'dap_input_manifest is required'
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    WRITE_COHORT_SUMMARY(file(params.dap_input_manifest, checkIfExists: true))
}
