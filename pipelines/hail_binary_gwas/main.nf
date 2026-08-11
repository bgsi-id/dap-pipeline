nextflow.enable.dsl = 2

include { HAIL_GWAS_CORE } from '../hail_gwas/main'

params.phenotype_type = 'binary'

workflow {
    HAIL_GWAS_CORE()
}
