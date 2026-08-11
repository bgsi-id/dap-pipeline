nextflow.enable.dsl = 2

include { HAIL_GWAS_CORE } from './main'

params.phenotype_type = 'quantitative'

workflow {
    HAIL_GWAS_CORE()
}
