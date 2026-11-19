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
// 'cohort': annotate every unannotated site in the variant table (default).
// 'samples': annotate only sites called in this run's samples (single-sample interpretation).
params.annotation_scope = 'cohort'
params.af_threshold = 0.01
params.indexed_annotate_max_sites = 50000
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
        sample_alias: (sample.source_alias ?: sample.display_name ?: sample.name ?: '')?.toString()?.trim(),
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
    require_file() {
      if [ ! -s "\$1" ]; then
        echo "Missing required annotation reference: \$1" >&2
        exit 1
      fi
    }
    require_index() {
      if [ ! -s "\$1.tbi" ] && [ ! -s "\$1.csi" ]; then
        echo "Missing tabix/CSI index for annotation reference: \$1" >&2
        exit 1
      fi
    }

    require_file '${params.reference_dir}/${params.fasta_name}'
    require_file '${params.reference_dir}/${params.fasta_name}.fai'
    require_file '${params.reference_dir}/${params.gff_name}'
    require_file '${params.reference_dir}/${params.gnomad_name}'
    require_file '${params.reference_dir}/${params.gnomad_popmax_name}'
    require_file '${params.reference_dir}/${params.clinvar_name}'
    require_index '${params.reference_dir}/${params.clinvar_name}'
    require_file '${params.reference_dir}/${params.spliceai_name}'
    require_index '${params.reference_dir}/${params.spliceai_name}'
    require_file '${params.reference_dir}/${params.revel_name}'
    require_file '${params.reference_dir}/${params.revel_name}.tbi'
    if [ ! -d '${params.vep_cache_dir}/homo_sapiens_merged/${params.vep_cache_version}_GRCh38' ]; then
      echo 'Missing required VEP cache: ${params.vep_cache_dir}/homo_sapiens_merged/${params.vep_cache_version}_GRCh38' >&2
      exit 1
    fi
    printf 'annotation_pack\t%s\nassembly\t%s\n' '${params.annotation_pack}' '${params.assembly}' > reference.ready
    """
}


process INGEST_SAMPLE {
    tag "${sample_id}"
    container params.variant_image
    cpus 4
    memory '12 GB'
    time '12h'
    publishDir "${params.dap_output_uri}/ingestion",
        mode: 'copy',
        overwrite: true,
        failOnError: true,
        pattern: 'ingested/*',
        saveAs: { name -> name.replaceFirst(/^ingested\//, '') }

    input:
    tuple val(sample_id), path(input_vcf), path(input_index)

    output:
    path "${sample_id}.meta.tsv", emit: meta
    path 'ingested/*', type: 'dir', emit: parquet_dir

    script:
    """
    set -euo pipefail
    vcf_sha256=\$(sha256sum '${input_vcf}' | awk '{print \$1}')
    run_id="${params.batch_id}-${sample_id}-\${vcf_sha256:0:16}"
    printf '%s\\t%s\\t%s\\n' '${sample_id}' "\${run_id}" "\${vcf_sha256}" > '${sample_id}.meta.tsv'
    mkdir -p "ingested/\${run_id}"

    python -m variant_ingest \
      --input '${input_vcf}' \
      --output "ingested/\${run_id}" \
      --release '${params.release_id}' \
      --batch '${params.batch_id}' \
      --assembly '${params.assembly}' \
      --sample-id '${sample_id}' \
      --allow-sample-id-mismatch \
      --run-id "\${run_id}" \
      --rows-per-file ${params.rows_per_file}
    """
}


process LOAD_VARIANTS {
    tag "${params.batch_id}"
    container params.variant_image
    cpus 1
    memory '2 GB'
    time '6h'
    maxForks 1
    errorStrategy 'retry'
    maxRetries 3

    input:
    path meta_files
    path parquet_trees

    output:
    path 'load-receipts.jsonl', emit: receipt

    script:
    """
    set -euo pipefail
    tab=\$(printf '\\t')
    : > load-receipts.jsonl
    for meta in *.meta.tsv; do
      IFS="\${tab}" read -r sample_id run_id vcf_sha256 < "\${meta}"
      variant_etl_control.py load-variants \
        --local-root "\${run_id}" \
        --published-root "${params.dap_output_uri}/ingestion/\${run_id}" \
        --clickhouse-url '${params.clickhouse_url}' \
        --database '${params.clickhouse_database}' \
        --region '${params.aws_region}' \
        --sample-id "\${sample_id}" \
        --source-sha256 "\${vcf_sha256}" \
        --output "\${sample_id}.load.json"
      cat "\${sample_id}.load.json" >> load-receipts.jsonl
      printf '\\n' >> load-receipts.jsonl
    done
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
    val scope_sample_ids

    output:
    path 'site.vcf.gz', emit: vcf
    path 'site-export.txt', emit: metrics

    script:
    def limitArg = params.max_annotation_variants ? "--max-variants ${params.max_annotation_variants}" : ''
    def scopeArg = scope_sample_ids ? "--sample-ids '${scope_sample_ids}'" : ''
    """
    set -euo pipefail
    CLICKHOUSE_URL='${params.clickhouse_url}' \
    CLICKHOUSE_DATABASE='${params.clickhouse_database}' \
    CLICKHOUSE_ANNOTATION_TABLE='${params.annotation_completion_table}' \
    python -m variant_ingest annotation-sites \
      --output site.vcf.gz \
      --annotation-pack '${params.annotation_pack}' \
      --assembly '${params.assembly}' \
      ${scopeArg} \
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
    path 'site.gnomad.raw.vcf.gz', emit: vcf

    script:
    """
    set -euo pipefail
    echtvar anno \
      -e '${params.reference_dir}/${params.gnomad_name}' \
      -e '${params.reference_dir}/${params.gnomad_popmax_name}' \
      '${sites_vcf}' site.gnomad.raw.vcf.gz
    """
}


process ANNOTATE_BASE_AND_SELECT {
    tag "${params.annotation_pack}"
    container params.bcftools_image
    cpus 8
    memory '16 GB'
    time '12h'

    input:
    path raw_vcf
    path reference_ready

    output:
    path 'site.base.vcf.gz', emit: base_vcf
    path 'site.base.vcf.gz.tbi', emit: base_index
    path 'site.small.vcf.gz', emit: small_vcf
    path 'site.small.vcf.gz.tbi', emit: small_index
    path 'selection.tsv', emit: metrics

    script:
    """
    set -euo pipefail
    tmp_dir=\$(mktemp -d)
    trap 'rm -rf "\${tmp_dir}"' EXIT

    header="\$(bcftools view -h '${raw_vcf}')"
    printf '%s\\n' "\${header}" | grep -q '##INFO=<ID=gnomad_af_max[>,]'
    : > "\${tmp_dir}/rename-annots.tsv"
    if ! printf '%s\\n' "\${header}" | grep -q '##INFO=<ID=gnomad_af_popmax[>,]'; then
      printf '%s\\n' "\${header}" | grep -q '##INFO=<ID=gnomad_popmax_af[>,]' \
        || { echo 'Echtvar output lacks gnomAD AF popmax field' >&2; exit 1; }
      printf 'INFO/gnomad_popmax_af INFO/gnomad_af_popmax\\n' >> "\${tmp_dir}/rename-annots.tsv"
    fi
    if ! printf '%s\\n' "\${header}" | grep -q '##INFO=<ID=gnomad_nhomalt[>,]'; then
      if printf '%s\\n' "\${header}" | grep -q '##INFO=<ID=gnomad_nhomalts[>,]'; then
        printf 'INFO/gnomad_nhomalts INFO/gnomad_nhomalt\\n' >> "\${tmp_dir}/rename-annots.tsv"
      fi
    fi
    if [ -s "\${tmp_dir}/rename-annots.tsv" ]; then
      bcftools annotate --rename-annots "\${tmp_dir}/rename-annots.tsv" \
        --threads ${task.cpus} -Oz -o "\${tmp_dir}/norm.vcf.gz" '${raw_vcf}'
    else
      cp '${raw_vcf}' "\${tmp_dir}/norm.vcf.gz"
    fi
    normalized_header="\$(bcftools view -h "\${tmp_dir}/norm.vcf.gz")"
    printf '%s\\n' "\${normalized_header}" | grep -q '##INFO=<ID=gnomad_af_popmax[>,]'
    bcftools index -t --threads ${task.cpus} "\${tmp_dir}/norm.vcf.gz"

    total=\$(bcftools index -n "\${tmp_dir}/norm.vcf.gz")
    region_args=""
    if [ "\${total}" -le ${params.indexed_annotate_max_sites} ]; then
      echo "Using indexed -R lookup for \${total} sites (threshold=${params.indexed_annotate_max_sites})" >&2
      region_args="-R \${tmp_dir}/norm.vcf.gz"
    fi

    bcftools annotate \${region_args} \
      -a '${params.reference_dir}/${params.clinvar_name}' \
      -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVI \
      --pair-logic exact \
      --threads ${task.cpus} \
      -Oz -o "\${tmp_dir}/clinvar.vcf.gz" \
      "\${tmp_dir}/norm.vcf.gz"
    bcftools index -t --threads ${task.cpus} "\${tmp_dir}/clinvar.vcf.gz"

    bcftools annotate \${region_args} \
      -a '${params.reference_dir}/${params.spliceai_name}' \
      -c INFO/SpliceAI \
      --pair-logic exact \
      -Ou "\${tmp_dir}/clinvar.vcf.gz" \
    | bcftools csq \
      --local-csq \
      --fasta-ref '${params.reference_dir}/${params.fasta_name}' \
      --gff-annot '${params.reference_dir}/${params.gff_name}' \
      --threads ${task.cpus} \
      -Oz -o site.base.vcf.gz
    bcftools index -t --threads ${task.cpus} site.base.vcf.gz

    bcftools view -Ov site.base.vcf.gz \
      | perl '${projectDir}/bin/add_rescue_fields.pl' \
      | bcftools view \
        --include '((INFO/gnomad_af_max="." || INFO/gnomad_af_max=-1 || INFO/gnomad_af_max<${params.af_threshold}) && (INFO/gnomad_af_popmax="." || INFO/gnomad_af_popmax=-1 || INFO/gnomad_af_popmax<${params.af_threshold})) || INFO/DAP_CLINVAR_RESCUE=1 || INFO/DAP_SPLICEAI_DS_MAX>=0.5' \
        --threads ${task.cpus} \
        -Oz -o site.small.vcf.gz
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
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple path(selected_vcf), path(selected_index), path(selection_metrics)
    path reference_ready

    output:
    tuple path('site.detail.vcf.gz'), path('site.detail.vcf.gz.tbi'), emit: detail

    script:
    """
    set -euo pipefail
    plugin_dir=""
    for candidate in "${projectDir}/bin" "\${VEP_DIR_PLUGINS:-/plugins}" /plugins /opt/vep/.vep/Plugins "${params.reference_dir}"; do
      if [ -f "\${candidate}/REVEL.pm" ]; then
        plugin_dir="\${candidate}"
        break
      fi
    done
    if [ -z "\${plugin_dir}" ]; then
      echo "REVEL.pm plugin not found in candidate paths" >&2
      exit 1
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
      --dir_plugins "\${plugin_dir}" \
      --plugin REVEL,file='${params.reference_dir}/${params.revel_name}' \
      --fork ${task.cpus} --buffer_size ${params.vep_buffer_size} \
      --no_stats --force_overwrite
    { zcat site.detail.vcf.gz || true; } | grep -m 1 '^##INFO=<ID=CSQ' | grep -qi 'REVEL'
    tabix -p vcf site.detail.vcf.gz
    """
}


process LOAD_ANNOTATIONS {
    tag "${params.annotation_pack}"
    container params.variant_image
    cpus 2
    memory '4 GB'
    time '4h'
    errorStrategy 'retry'
    maxRetries 2
    publishDir params.dap_output_uri, mode: 'copy', overwrite: true, failOnError: true

    input:
    path 'site.base.vcf.gz'
    path 'site.base.vcf.gz.tbi'
    tuple path('site.detail.vcf.gz'), path('site.detail.vcf.gz.tbi')
    path selection_metrics
    path site_export_metrics

    output:
    path 'annotation-load.json', emit: receipt
    path 'results', emit: published

    script:
    """
    set -euo pipefail
    variant_etl_control.py load-annotations \
      --base-vcf site.base.vcf.gz \
      --detail-vcf site.detail.vcf.gz \
      --annotation-pack '${params.annotation_pack}' \
      --assembly '${params.assembly}' \
      --clickhouse-url '${params.clickhouse_url}' \
      --database '${params.clickhouse_database}' \
      --completion-table '${params.annotation_completion_table}' \
      --output annotation-load.json

    mkdir -p results/annotation
    cp site.base.vcf.gz site.base.vcf.gz.tbi results/annotation/
    cp site.detail.vcf.gz site.detail.vcf.gz.tbi results/annotation/
    cp '${selection_metrics}' annotation-load.json '${site_export_metrics}' results/annotation/
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
    if (!(params.annotation_scope in ['cohort', 'samples'])) {
        error "annotation_scope must be 'cohort' or 'samples'"
    }
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
    references = PRECHECK_REFERENCES()
    ingested = INGEST_SAMPLE(sample_inputs)
    loaded = LOAD_VARIANTS(ingested.meta.collect(), ingested.parquet_dir.collect())
    def scopeIds = (
        sampleIds.collect { safe_component(it, 'sample_id') } +
        records.collect { it.sample_alias }.findAll { it && (it ==~ /[A-Za-z0-9._-]{1,128}/) }
    ).unique()
    scope_sample_ids = params.annotation_scope == 'samples' ? scopeIds.join(',') : ''
    sites = EXPORT_NOVEL_SITES(loaded.receipt, references, scope_sample_ids)
    af_raw = ANNOTATE_AF(sites.vcf, references)
    base = ANNOTATE_BASE_AND_SELECT(af_raw.vcf, references)
    detail_branch = base.small_vcf
        .combine(base.small_index)
        .combine(base.metrics)
        .branch { vcf, idx, metrics ->
            def selectedCount = metrics.readLines().findResult { line ->
                def parts = line.split('\t')
                parts.length >= 2 && parts[0] == 'selected_sites' ? parts[1].trim().toLong() : null
            }
            if (selectedCount == null) error 'selection.tsv missing selected_sites'
            run_vep: selectedCount > 0
                return tuple(vcf, idx, metrics)
            skip_vep: true
                return tuple(vcf, idx)
        }
    vep = ANNOTATE_DETAIL_VEP(detail_branch.run_vep, references)
    detail_out = vep.detail.mix(detail_branch.skip_vep)
    LOAD_ANNOTATIONS(base.base_vcf, base.base_index, detail_out, base.metrics, sites.metrics)
}
