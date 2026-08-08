# DAP Nextflow foundation

`base_nf` is the small contract layer used by DAP downstream workflows. It is
not an nf-core replacement and it does not resolve authorization or inspect
governed data itself.

`dap-runtime` supplies an immutable `urn:bgsi:dap:resolved-inputs:2` document.
The document contains named cohort memberships and release-level asset URIs
that have already been authorized by DAP Data.

## Modules

`READ_BUNDLES` accepts the resolved manifest as a `path` channel and emits one
normalized `bundle` value:

```nextflow
include { READ_BUNDLES } from '../base_nf/modules/bundle'

bundles = READ_BUNDLES(Channel.value(file(params.dap_input_manifest)))
bundles.bundle.map { bundle ->
    def control = bundle.cohorts.control.samples
    def pvcf = bundle.release_assets.pvcf.access_uri
}
```

The module validates the manifest schema, dataset/release identity, cohort
sample lists, and release-asset roles. It does not decide which cohort names or
asset roles a scientific workflow needs.

`WRITE_RESULTS` accepts a final result directory and copies it to the project
output prefix in `params.dap_output_uri`. Publishing failures fail the
workflow.

```nextflow
include { WRITE_RESULTS } from '../base_nf/modules/bundle'

WRITE_RESULTS(final_results)
```

## Foundation smoke workflow

`main.nf` validates the DAP contract and writes a non-sensitive bundle receipt
with only dataset identity, cohort counts, and asset roles. It is useful for a
WES integration check; it does not process release objects.

Every downstream pipeline should keep its own scientific steps explicit and
small. It should use `READ_BUNDLES` at the boundary and `WRITE_RESULTS` once
for declared final outputs.
