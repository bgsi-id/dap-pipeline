nextflow.enable.dsl = 2

/*
 * Cohort variant ETL
 *
 * governed VCFs -> Parquet facts -> ClickHouse variant dimension
 *   -> never-annotated sites -> gnomAD AF -> local bcftools CSQ
 *   -> rare/missing-AF sites -> VEP -> ClickHouse annotation projections
 *
 * All object-store inputs and published outputs are handled by Nextflow.
 * Processes do not invoke an object-store CLI.
 */

params.dap_input_manifest = null
params.dap_output_uri = null
params.output_dir = 'results'
params.release_id = null
params.batch_id = null
params.assembly = 'GRCh38'
params.annotation_pack = 'grch38-v3'
params.af_threshold = 0.01
params.max_annotation_variants = null
params.reference_dir = '/reference'
params.fasta_name = 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
params.gff_name = 'Homo_sapiens.GRCh38.116.chr.gff3.gz'
params.gnomad_name = 'gnomad_v4.1.zip'
params.gnomad_popmax_name = 'gnomad.v3.1.2.echtvar.popmax.v2.zip'
params.clinvar_name = 'clinvar.chr.vcf.gz'
params.spliceai_name = 'spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz'
params.revel_name = 'revel_grch38.tsv.gz'
params.vep_cache_dir = '/reference/vep'
params.vep_cache_version = '116'
params.vep_buffer_size = 5000
params.rows_per_file = 50000
params.aws_region = 'ap-southeast-3'
params.clickhouse_url = 'http://dap-ch:8123'
params.clickhouse_database = 'default'
params.annotation_completion_table = 'variant_annotation_complete'


def safe_component(def value, String label) {
    def text = value?.toString()?.trim()
    if (!text || !(text ==~ /[A-Za-z0-9._-]{1,128}/)) {
        error "${label} must contain only A-Z, a-z, 0-9, dot, underscore or hyphen"
    }
    return text
}


def sample_from_asset(def sample) {
    def assets = (sample.assets ?: []).collectEntries { asset ->
        [(asset.role?.toString()?.toLowerCase()): asset]
    }
    if (!sample.id || !assets.vcf || !assets.vcf_index) {
        error 'Each resolved sample requires id, vcf and vcf_index assets'
    }
    return [
        sample_id: sample.id.toString(),
        vcf_uri: assets.vcf.access_uri.toString(),
        index_uri: assets.vcf_index.access_uri.toString(),
    ]
}


process PRECHECK_REFERENCES {
    tag "${params.annotation_pack}"
    container params.bcftools_image
    cpus 1
    memory '2 GB'
    time '30m'

    output:
    path 'reference.ready'

    script:
    """
    set -euo pipefail
    test -s '${params.reference_dir}/${params.fasta_name}'
    test -s '${params.reference_dir}/${params.fasta_name}.fai'
    test -s '${params.reference_dir}/${params.gff_name}'
    test -s '${params.reference_dir}/${params.gnomad_name}'
    test -s '${params.reference_dir}/${params.gnomad_popmax_name}'
    test -s '${params.reference_dir}/${params.clinvar_name}'
    test -s '${params.reference_dir}/${params.clinvar_name}.tbi' -o -s '${params.reference_dir}/${params.clinvar_name}.csi'
    test -s '${params.reference_dir}/${params.spliceai_name}'
    test -s '${params.reference_dir}/${params.spliceai_name}.tbi' -o -s '${params.reference_dir}/${params.spliceai_name}.csi'
    test -s '${params.reference_dir}/${params.revel_name}'
    test -s '${params.reference_dir}/${params.revel_name}.tbi'
    test -d '${params.vep_cache_dir}/homo_sapiens_merged/${params.vep_cache_version}_GRCh38'
    printf 'annotation_pack\t%s\nassembly\t%s\n' '${params.annotation_pack}' '${params.assembly}' > reference.ready
    """
}


process PREPARE_SAMPLE_INPUT {
    tag "${sample_id}"
    container params.bcftools_image
    cpus 1
    memory '2 GB'
    time '2h'

    input:
    tuple val(sample_id), path(input_vcf), path(input_index)

    output:
    tuple val(sample_id), path('input.vcf.gz'), path('input.vcf.gz.tbi'), path('checksums.json')

    script:
    """
    set -euo pipefail
    cp '${input_vcf}' input.vcf.gz
    cp '${input_index}' input.vcf.gz.tbi
    vcf_sha256=\$(sha256sum input.vcf.gz | awk '{print \$1}')
    index_sha256=\$(sha256sum input.vcf.gz.tbi | awk '{print \$1}')
    printf '{"vcf_sha256":"%s","index_sha256":"%s"}\\n' "\${vcf_sha256}" "\${index_sha256}" > checksums.json
    """
}


process INGEST_SAMPLE {
    tag "${sample_id}"
    container params.variant_image
    cpus 4
    memory '12 GB'
    time '12h'
    publishDir "${params.dap_output_uri}/ingestion", mode: 'copy', overwrite: true, failOnError: true

    input:
    tuple val(sample_id), path(input_vcf), val(vcf_sha256), path(input_index), val(index_sha256), val(run_id)

    output:
    tuple val(sample_id), val(run_id), val(vcf_sha256), path("${run_id}"), emit: parquet

    script:
    """
    set -euo pipefail
    echo '${vcf_sha256}  ${input_vcf}' | sha256sum -c -
    echo '${index_sha256}  ${input_index}' | sha256sum -c -

    python -m variant_ingest \
      --input '${input_vcf}' \
      --output '${run_id}' \
      --release '${params.release_id}' \
      --batch '${params.batch_id}' \
      --assembly '${params.assembly}' \
      --sample-id '${sample_id}' \
      --allow-sample-id-mismatch \
      --run-id '${run_id}' \
      --rows-per-file ${params.rows_per_file}
    """
}


process LOAD_VARIANTS {
    tag "${sample_id}"
    container params.variant_image
    cpus 1
    memory '2 GB'
    time '2h'
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample_id), val(run_id), val(vcf_sha256), path(parquet_tree)

    output:
    tuple val(sample_id), path("${sample_id}.load.json"), emit: receipt

    script:
    """
    set -euo pipefail
    variant_etl_control.py load-variants \
      --local-root '${parquet_tree}' \
      --published-root '${params.dap_output_uri}/ingestion/${run_id}' \
      --clickhouse-url '${params.clickhouse_url}' \
      --database '${params.clickhouse_database}' \
      --region '${params.aws_region}' \
      --sample-id '${sample_id}' \
      --source-sha256 '${vcf_sha256}' \
      --output '${sample_id}.load.json'
    """
}


process EXPORT_NOVEL_SITES {
    tag "${params.annotation_pack}"
    container params.variant_image
    cpus 2
    memory '4 GB'
    time '4h'

    input:
    val load_receipts
    path reference_ready

    output:
    path 'site.vcf.gz', emit: vcf
    path 'site-export.txt', emit: metrics

    script:
    def limitArg = params.max_annotation_variants ? "--max-variants ${params.max_annotation_variants}" : ''
    """
    set -euo pipefail
    CLICKHOUSE_URL='${params.clickhouse_url}' \
    CLICKHOUSE_DATABASE='${params.clickhouse_database}' \
    CLICKHOUSE_ANNOTATION_TABLE='${params.annotation_completion_table}' \
    python -m variant_ingest annotation-sites \
      --output site.vcf.gz \
      --annotation-pack '${params.annotation_pack}' \
      --assembly '${params.assembly}' \
      ${limitArg} | tee site-export.txt
    """
}


process ANNOTATE_AF {
    tag "${params.annotation_pack}"
    container params.echtvar_image
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    path sites_vcf
    path reference_ready

    output:
    path 'site.gnomad.vcf.gz', emit: vcf

    script:
    """
    set -euo pipefail
    echtvar anno \
      -e '${params.reference_dir}/${params.gnomad_name}' \
      -e '${params.reference_dir}/${params.gnomad_popmax_name}' \
      '${sites_vcf}' site.gnomad.vcf.gz
    """
}


process ANNOTATE_CLINVAR {
    tag "${params.annotation_pack}"
    container params.bcftools_image
    cpus 8
    memory '12 GB'
    time '4h'

    input:
    path frequency_vcf
    path reference_ready

    output:
    path 'site.clinvar.vcf.gz', emit: vcf
    path 'site.clinvar.vcf.gz.tbi', emit: index

    script:
    """
    set -euo pipefail
    bcftools index -t --threads ${task.cpus} '${frequency_vcf}'
    bcftools annotate \
      -a '${params.reference_dir}/${params.clinvar_name}' \
      -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVI \
      --pair-logic exact \
      --threads ${task.cpus} \
      -Oz -o site.clinvar.vcf.gz \
      '${frequency_vcf}'
    bcftools index -t --threads ${task.cpus} site.clinvar.vcf.gz
    """
}


process ANNOTATE_SPLICEAI {
    tag "${params.annotation_pack}"
    container params.bcftools_image
    cpus 8
    memory '12 GB'
    time '4h'

    input:
    path clinvar_vcf
    path clinvar_index
    path reference_ready

    output:
    path 'site.spliceai.vcf.gz', emit: vcf
    path 'site.spliceai.vcf.gz.tbi', emit: index

    script:
    """
    set -euo pipefail
    bcftools annotate \
      -a '${params.reference_dir}/${params.spliceai_name}' \
      -c INFO/SpliceAI \
      --pair-logic exact \
      --threads ${task.cpus} \
      -Oz -o site.spliceai.vcf.gz \
      '${clinvar_vcf}'
    bcftools index -t --threads ${task.cpus} site.spliceai.vcf.gz
    """
}


process ANNOTATE_LOCAL_CSQ {
    tag "${params.annotation_pack}"
    container params.bcftools_image
    cpus 8
    memory '16 GB'
    time '8h'

    input:
    path annotated_vcf
    path annotated_index
    path reference_ready

    output:
    path 'site.base.vcf.gz', emit: vcf
    path 'site.base.vcf.gz.tbi', emit: index

    script:
    """
    set -euo pipefail
    bcftools csq \
      --local-csq \
      --fasta-ref '${params.reference_dir}/${params.fasta_name}' \
      --gff-annot '${params.reference_dir}/${params.gff_name}' \
      --threads ${task.cpus} \
      -Oz -o site.base.vcf.gz \
      '${annotated_vcf}'
    bcftools index -t --threads ${task.cpus} site.base.vcf.gz
    """
}


process SELECT_DETAIL_SITES {
    tag "AF<${params.af_threshold}"
    container params.bcftools_image
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    path base_vcf
    path base_index

    output:
    path 'site.small.vcf.gz', emit: vcf
    path 'site.small.vcf.gz.tbi', emit: index
    path 'selection.tsv', emit: metrics

    script:
    """
    set -euo pipefail
    total=\$(bcftools index -n '${base_vcf}')
    bcftools view \
      --include '((INFO/gnomad_af_max="." || INFO/gnomad_af_max=-1 || INFO/gnomad_af_max<${params.af_threshold}) && (INFO/gnomad_af_popmax="." || INFO/gnomad_af_popmax=-1 || INFO/gnomad_af_popmax<${params.af_threshold})) || INFO/CLNSIG ~ "(?i)(pathogenic|likely_pathogenic)" || INFO/CLNSIGCONF ~ "(?i)(pathogenic|likely_pathogenic)" || INFO/SpliceAI ~ "(0\\.[5-9]|1\\.0)"' \
      --threads ${task.cpus} \
      -Oz -o site.small.vcf.gz \
      '${base_vcf}'
    bcftools index -t --threads ${task.cpus} site.small.vcf.gz
    selected=\$(bcftools index -n site.small.vcf.gz)
    printf 'total_sites\\t%s\\nselected_sites\\t%s\\naf_threshold\\t%s\\n' \
      "\${total}" "\${selected}" '${params.af_threshold}' > selection.tsv
    """
}


process ANNOTATE_DETAIL_VEP {
    tag "${params.annotation_pack}"
    container params.vep_image
    cpus 16
    memory '32 GB'
    time '18h'

    input:
    path selected_vcf
    path selected_index
    path selection_metrics
    path reference_ready

    output:
    path 'site.detail.vcf.gz', emit: vcf

    script:
    """
    set -euo pipefail
    selected_count=\$(sed -n 's/^selected_sites[[:space:]]*//p' '${selection_metrics}')
    test -n "\${selected_count}"
    if [ "\${selected_count}" -eq 0 ]; then
      cp '${selected_vcf}' site.detail.vcf.gz
    else
      revel_arg=""
      if [ -f '${params.reference_dir}/${params.revel_name}' ]; then
        revel_arg="--plugin REVEL,file=${params.reference_dir}/${params.revel_name}"
      fi
      vep \
        --input_file '${selected_vcf}' \
        --output_file site.detail.vcf.gz \
        --vcf --compress_output bgzip --offline --cache \
        --dir_cache '${params.vep_cache_dir}' \
        --merged --cache_version '${params.vep_cache_version}' \
        --assembly '${params.assembly}' \
        --fasta '${params.reference_dir}/${params.fasta_name}' \
        --mane --canonical --symbol --biotype --hgvs --hgvsg \
        --allele_number --sift b --polyphen b \
        --shift_hgvs 1 --numbers --domains --protein --uniprot \
        --flag_pick --pick_order mane_select,mane_plus_clinical,canonical,rank \
        \${revel_arg} \
        --fork ${task.cpus} --buffer_size ${params.vep_buffer_size} \
        --no_stats --force_overwrite
    fi
    """
}


process INDEX_DETAIL_VCF {
    tag "${params.annotation_pack}"
    container params.bcftools_image
    cpus 4
    memory '4 GB'
    time '2h'

    input:
    path vep_vcf, name: 'vep-output.vcf.gz'

    output:
    path 'site.detail.vcf.gz', emit: vcf
    path 'site.detail.vcf.gz.tbi', emit: index

    script:
    """
    set -euo pipefail
    cp vep-output.vcf.gz site.detail.vcf.gz
    bcftools index -t --threads ${task.cpus} site.detail.vcf.gz
    """
}


process LOAD_ANNOTATIONS {
    tag "${params.annotation_pack}"
    container params.variant_image
    cpus 2
    memory '4 GB'
    time '4h'

    input:
    path base_vcf
    path detail_vcf

    output:
    path 'annotation-load.json', emit: receipt

    script:
    """
    set -euo pipefail
    variant_etl_control.py load-annotations \
      --base-vcf '${base_vcf}' \
      --detail-vcf '${detail_vcf}' \
      --annotation-pack '${params.annotation_pack}' \
      --assembly '${params.assembly}' \
      --clickhouse-url '${params.clickhouse_url}' \
      --database '${params.clickhouse_database}' \
      --completion-table '${params.annotation_completion_table}' \
      --output annotation-load.json
    """
}


process COLLECT_RESULTS {
    tag "${params.batch_id}"
    container params.variant_image
    cpus 1
    memory '2 GB'
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true, failOnError: true

    input:
    path site_base_vcf
    path site_base_index
    path site_detail_vcf
    path site_detail_index
    path selection_metrics
    path annotation_receipt
    path site_export_metrics

    output:
    path 'results'

    script:
    """
    set -euo pipefail
    mkdir -p results/annotation
    cp '${site_base_vcf}' '${site_base_index}' results/annotation/
    cp '${site_detail_vcf}' '${site_detail_index}' results/annotation/
    cp '${selection_metrics}' '${annotation_receipt}' '${site_export_metrics}' results/annotation/
    python - <<'PY'
import json, pathlib
out = pathlib.Path('results/run.json')
out.write_text(json.dumps({
    'pipeline': 'variant-etl',
    'release_id': '${params.release_id}',
    'batch_id': '${params.batch_id}',
    'assembly': '${params.assembly}',
    'annotation_pack': '${params.annotation_pack}',
    'af_threshold': ${params.af_threshold},
}, indent=2, sort_keys=True))
PY
    """
}


workflow {
    if (!params.dap_input_manifest) error 'dap_input_manifest is required'
    if (!params.dap_output_uri) error 'dap_output_uri is required'
    if (!params.variant_image) {
        error 'variant_image is required and must be injected by the deployment/runtime'
    }
    safe_component(params.release_id, 'release_id')
    safe_component(params.batch_id, 'batch_id')
    safe_component(params.annotation_pack, 'annotation_pack')
    def threshold = params.af_threshold as double
    if (threshold < 0.0 || threshold > 1.0) error 'af_threshold must be between 0 and 1'

    def resolved = new groovy.json.JsonSlurper().parse(file(params.dap_input_manifest).toFile())
    if (resolved.schema != 'urn:bgsi:dap:resolved-inputs:1' || !(resolved.samples instanceof List)) {
        error 'dap_input_manifest must use resolved-inputs v1 with a samples list'
    }
    def records = resolved.samples.collect { sample_from_asset(it) }
    if (!records) error 'No input samples were resolved'

    def sampleIds = records.collect { (it.sample_id ?: it.id)?.toString() }
    if (sampleIds.any { !it } || sampleIds.toSet().size() != sampleIds.size()) {
        error 'Input sample IDs must be present and unique'
    }

    sample_inputs = Channel.fromList(records).map { record ->
        def sampleId = safe_component(record.sample_id ?: record.id, 'sample_id')
        def vcfUri = record.vcf_uri?.toString()
        def indexUri = record.index_uri?.toString()
        if (!vcfUri || !indexUri) {
            error "Sample ${sampleId} requires VCF and VCF index URIs"
        }
        tuple(sampleId, file(vcfUri, checkIfExists: true), file(indexUri, checkIfExists: true))
    }
    prepared_inputs = PREPARE_SAMPLE_INPUT(sample_inputs).map { sampleId, vcf, index, checksumFile ->
        def checksums = new groovy.json.JsonSlurper().parse(checksumFile.toFile())
        def sha = checksums.vcf_sha256?.toString()?.toLowerCase()
        def indexSha = checksums.index_sha256?.toString()?.toLowerCase()
        if (!(sha ==~ /[0-9a-f]{64}/) || !(indexSha ==~ /[0-9a-f]{64}/)) {
            error "Sample ${sampleId} checksum preparation failed"
        }
        def runId = safe_component("${params.batch_id}-${sampleId}-${sha.take(16)}", 'run_id')
        tuple(sampleId, vcf, sha, index, indexSha, runId)
    }

    references = PRECHECK_REFERENCES()
    ingested = INGEST_SAMPLE(prepared_inputs)
    loaded = LOAD_VARIANTS(ingested.parquet)
    load_barrier = loaded.receipt.collect()
    sites = EXPORT_NOVEL_SITES(load_barrier, references)
    af = ANNOTATE_AF(sites.vcf, references)
    clinvar = ANNOTATE_CLINVAR(af.vcf, references)
    spliceai = ANNOTATE_SPLICEAI(clinvar.vcf, clinvar.index, references)
    base = ANNOTATE_LOCAL_CSQ(spliceai.vcf, spliceai.index, references)
    selected = SELECT_DETAIL_SITES(base.vcf, base.index)
    vep = ANNOTATE_DETAIL_VEP(selected.vcf, selected.index, selected.metrics, references)
    detail = INDEX_DETAIL_VCF(vep.vcf)
    annotation_load = LOAD_ANNOTATIONS(base.vcf, detail.vcf)
    COLLECT_RESULTS(base.vcf, base.index, detail.vcf, detail.index,
                    selected.metrics, annotation_load.receipt, sites.metrics)
}
