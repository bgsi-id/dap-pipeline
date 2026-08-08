nextflow.enable.dsl = 2

/*
 * Rare-disease WGS clinical annotation
 *
 * normalize -> gnomAD -> ClinVar -> prefilter -> VEP -> review TSV
 *
 * Reference data is mounted read-only at a stable path by the execution
 * environment. It is not staged or copied by Nextflow.
 */

params.vcf_uri = null
params.vcf_index_uri = null
params.vcf_sha256 = ''
params.vcf_index_sha256 = ''
params.sample_id = null
params.reference_release = null
params.reference_dir = '/reference'
params.output_dir = 'results'
params.dap_output_uri = null
params.assay_type = 'wgs'
params.af_cutoff = 0.01
params.af_max_cutoff = 0.02

params.fasta_name = 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
params.gnomad_name = 'gnomad_v4.1.zip'
params.clinvar_name = 'clinvar.chr.vcf.gz'
params.spliceai_name = 'spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz'
params.spliceai_indel_name = ''
params.vep_cache_subdir = 'auto'
params.vep_cache_version = '116'
params.vep_buffer_size = 5000

params.bcftools_image = 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
params.echtvar_image = 'quay.io/biocontainers/echtvar@sha256:71cd0028f4aa9f7d012be3dc86b20e4df78fc5d8178d323c3c215a3f449a6244'
params.vep_image = 'ensemblorg/ensembl-vep:release_116.0'


process PRECHECK_REFERENCE {
    tag "${reference_release}"

    container params.bcftools_image
    cpus 1
    memory '1 GB'
    time '30m'

    input:
    val reference_release

    output:
    path 'reference.ready', emit: ready

    script:
    """
    set -euo pipefail

    reference_dir='${params.reference_dir}'
    fasta="\${reference_dir}/${params.fasta_name}"
    clinvar="\${reference_dir}/${params.clinvar_name}"

    test -s "\${fasta}"
    test -s "\${fasta}.fai"
    test -s "\${reference_dir}/${params.gnomad_name}"
    test -s "\${clinvar}"
    test -s "\${clinvar}.tbi"

    cut -f1 "\${fasta}.fai" | LC_ALL=C sort -u > fasta-contigs.txt
    bcftools index -s "\${clinvar}" | cut -f1 | LC_ALL=C sort -u > clinvar-contigs.txt
    comm -12 fasta-contigs.txt clinvar-contigs.txt > shared-contigs.txt
    if [[ ! -s shared-contigs.txt ]]; then
      echo 'ERROR: ClinVar and FASTA have no contigs in common' >&2
      exit 1
    fi

    configured_cache='${params.vep_cache_subdir}'
    vep_cache_dir=''
    if [[ "\${configured_cache}" != 'auto' ]]; then
      candidate="\${reference_dir}/\${configured_cache}"
      if [[ -d "\${candidate}/homo_sapiens_merged/${params.vep_cache_version}_GRCh38" ]]; then
        vep_cache_dir="\${candidate}"
      fi
    else
      for candidate in \
        "\${reference_dir}/vep" \
        "\${reference_dir}/homo_sapiens_merged_vep_${params.vep_cache_version}_GRCh38" \
        "\${reference_dir}"
      do
        if [[ -d "\${candidate}/homo_sapiens_merged/${params.vep_cache_version}_GRCh38" ]]; then
          vep_cache_dir="\${candidate}"
          break
        fi
      done
    fi
    if [[ -z "\${vep_cache_dir}" ]]; then
      echo 'ERROR: VEP merged cache hierarchy homo_sapiens_merged/${params.vep_cache_version}_GRCh38 not found' >&2
      exit 1
    fi

    if [[ -n '${params.spliceai_indel_name}' ]]; then
      test -s "\${reference_dir}/${params.spliceai_name}"
      test -s "\${reference_dir}/${params.spliceai_name}.tbi"
      test -s "\${reference_dir}/${params.spliceai_indel_name}"
      test -s "\${reference_dir}/${params.spliceai_indel_name}.tbi"
    fi

    printf 'reference_release\\t%s\\nvep_cache_dir\\t%s\\n' \
      '${reference_release}' "\${vep_cache_dir}" > reference.ready
    """
}


process FILTER_REFERENCE_CONTIGS {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    path input_vcf
    path input_vcf_index
    path reference_ready
    val sample_id
    val vcf_sha256
    val vcf_index_sha256
    val reference_release

    output:
    path 'primary.vcf.gz', emit: vcf
    path 'primary.vcf.gz.tbi', emit: index
    path 'contig-filter.tsv', emit: stats

    script:
    """
    set -euo pipefail

    if [[ -n '${vcf_sha256}' ]]; then
      echo '${vcf_sha256}  ${input_vcf}' | sha256sum -c -
    fi
    if [[ -n '${vcf_index_sha256}' ]]; then
      echo '${vcf_index_sha256}  ${input_vcf_index}' | sha256sum -c -
    fi

    awk 'BEGIN { OFS="\\t" } { print \$1, 1, \$2 }' \
      '${params.reference_dir}/${params.fasta_name}.fai' \
      > reference-contigs.tsv

    bcftools view \
      --regions-file reference-contigs.tsv \
      --threads ${task.cpus} \
      -Oz \
      -o primary.vcf.gz \
      '${input_vcf}'
    bcftools index -t --threads ${task.cpus} primary.vcf.gz

    variants_total=\$(bcftools index -n '${input_vcf}')
    variants_retained=\$(bcftools index -n primary.vcf.gz)
    variants_removed=\$((variants_total - variants_retained))
    if (( variants_retained == 0 )); then
      echo 'ERROR: no variants remain after filtering to FASTA contigs' >&2
      exit 1
    fi
    printf 'variants_total\\t%s\\nvariants_retained\\t%s\\nvariants_removed_non_reference_contigs\\t%s\\n' \
      "\${variants_total}" "\${variants_retained}" "\${variants_removed}" \
      > contig-filter.tsv
    """
}


process NORMALIZE {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 16
    memory '16 GB'
    time '8h'

    input:
    path input_vcf
    path input_vcf_index
    path reference_ready
    val sample_id
    val reference_release

    output:
    path 'norm.vcf.gz', emit: vcf
    path 'norm.vcf.gz.tbi', emit: index

    script:
    """
    set -euo pipefail

    bcftools norm \
      -m -any \
      -f '${params.reference_dir}/${params.fasta_name}' \
      -Ou '${input_vcf}' \
    | bcftools norm \
      -d exact \
      --threads ${task.cpus} \
      -Oz \
      -o norm.vcf.gz

    bcftools index -t --threads ${task.cpus} norm.vcf.gz
    """
}


process ANNOTATE_GNOMAD {
    tag "${sample_id}"

    container params.echtvar_image
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    path normalized_vcf
    path normalized_index
    path reference_ready
    val sample_id
    val reference_release

    output:
    path 'freq.vcf.gz', emit: vcf

    script:
    """
    set -euo pipefail

    echtvar anno \
      -e '${params.reference_dir}/${params.gnomad_name}' \
      '${normalized_vcf}' \
      freq.vcf.gz
    """
}


process ANNOTATE_CLINVAR {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 8
    memory '12 GB'
    time '4h'

    input:
    path frequency_vcf
    path reference_ready
    val sample_id
    val reference_release

    output:
    path 'annot.vcf.gz', emit: vcf
    path 'annot.vcf.gz.tbi', emit: index
    path 'annotation-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    bcftools index -t --threads ${task.cpus} '${frequency_vcf}'
    bcftools annotate \
      -a '${params.reference_dir}/${params.clinvar_name}' \
      -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVI \
      --pair-logic exact \
      --threads ${task.cpus} \
      -Oz \
      -o annot.vcf.gz \
      '${frequency_vcf}'
    bcftools index -t --threads ${task.cpus} annot.vcf.gz

    bcftools query \
      -f '%INFO/gnomad_af\\t%INFO/gnomad_af_max\\t%INFO/CLNSIG\\t%INFO/CLNSIGCONF\\n' \
      annot.vcf.gz \
    | awk -F '\\t' '
        BEGIN { OFS="\\t" }
        {
          total++
          if (\$1 != "." && \$1 != "-1") gnomad_af_n++
          if (\$2 != "." && \$2 != "-1") gnomad_af_max_n++
          if (\$3 != ".") {
            clinvar_n++
            value=tolower(\$3)
            if (value ~ /(^|[,\\/|])(pathogenic|likely_pathogenic)([,\\/|]|\$)/) clinvar_plp_n++
            if (value ~ /conflicting_classifications_of_pathogenicity/) clinvar_conflict_n++
          }
          if (\$4 != ".") clinvar_conf_detail_n++
        }
        END {
          print "annotation_variants", total+0
          print "gnomad_af_annotated", gnomad_af_n+0
          print "gnomad_af_max_annotated", gnomad_af_max_n+0
          print "clinvar_annotated", clinvar_n+0
          print "clinvar_plp", clinvar_plp_n+0
          print "clinvar_conflicting", clinvar_conflict_n+0
          print "clinvar_conflict_detail", clinvar_conf_detail_n+0
        }
      ' > annotation-metrics.tsv

    clinvar_annotated=\$(awk -F '\\t' '\$1=="clinvar_annotated" { print \$2 }' annotation-metrics.tsv)
    if [[ '${params.assay_type}' == 'wgs' && "\${clinvar_annotated}" -eq 0 ]]; then
      echo 'ERROR: ClinVar annotation matched zero records for a WGS input' >&2
      exit 1
    fi
    """
}


process PREFILTER {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 8
    memory '12 GB'
    time '4h'

    input:
    path annotated_vcf
    path annotated_index
    path annotation_metrics
    val sample_id
    val af_cutoff
    val af_max_cutoff

    output:
    path 'tiered.vcf.gz', emit: vcf
    path 'tiered.vcf.gz.tbi', emit: index
    path 'variant-counts.tsv', emit: counts

    script:
    """
    set -euo pipefail

    bcftools filter \
      -i '((INFO/gnomad_af = "." || INFO/gnomad_af = -1 || INFO/gnomad_af < ${af_cutoff}) && (INFO/gnomad_af_max = "." || INFO/gnomad_af_max = -1 || INFO/gnomad_af_max < ${af_max_cutoff})) || INFO/CLNSIG ~ "athogenic" || INFO/CLNSIGCONF != "."' \
      --threads ${task.cpus} \
      -Oz \
      -o tiered.vcf.gz \
      '${annotated_vcf}'
    bcftools index -t --threads ${task.cpus} tiered.vcf.gz

    variants_in=\$(bcftools index -n '${annotated_vcf}')
    variants_tiered=\$(bcftools index -n tiered.vcf.gz)
    cp '${annotation_metrics}' variant-counts.tsv
    printf 'variants_in\\t%s\\nvariants_tiered\\t%s\\n' \
      "\${variants_in}" "\${variants_tiered}" >> variant-counts.tsv

    bcftools query -f '%INFO/CLNSIG\\t%INFO/CLNSIGCONF\\n' tiered.vcf.gz \
    | awk -F '\\t' '
        BEGIN { OFS="\\t" }
        {
          value=tolower(\$1)
          if (value ~ /(^|[,\\/|])(pathogenic|likely_pathogenic)([,\\/|]|\$)/) plp_n++
          if (value ~ /conflicting_classifications_of_pathogenicity/ || \$2 != ".") conflict_n++
        }
        END {
          print "prefilter_plp_retained", plp_n+0
          print "prefilter_conflicting_retained", conflict_n+0
        }
      ' >> variant-counts.tsv
    """
}


process ANNOTATE_VEP {
    tag "${sample_id}"

    container params.vep_image
    cpus 16
    memory '32 GB'
    time '12h'

    input:
    path tiered_vcf
    path tiered_index
    path reference_ready
    val sample_id
    val reference_release

    output:
    path "${sample_id}.annotated.vcf.gz", emit: vcf
    path 'vep-run-metadata.tsv', emit: metadata

    script:
    """
    set -euo pipefail

    vep_cache_dir=\$(awk -F '\\t' '\$1=="vep_cache_dir" { print \$2 }' '${reference_ready}')
    if [[ -z "\${vep_cache_dir}" ]]; then
      echo 'ERROR: resolved VEP cache directory is missing from reference.ready' >&2
      exit 1
    fi

    spliceai_args=()
    spliceai_status='NOT_APPLIED'
    if [[ -n '${params.spliceai_indel_name}' ]]; then
      spliceai_args+=(
        --plugin
        'SpliceAI,snv=${params.reference_dir}/${params.spliceai_name},indel=${params.reference_dir}/${params.spliceai_indel_name}'
      )
      spliceai_status='APPLIED'
    fi

    vep \
      --input_file '${tiered_vcf}' \
      --output_file '${sample_id}.annotated.vcf.gz' \
      --vcf \
      --compress_output bgzip \
      --offline \
      --cache \
      --dir_cache "\${vep_cache_dir}" \
      --merged \
      --cache_version '${params.vep_cache_version}' \
      --assembly GRCh38 \
      --fasta '${params.reference_dir}/${params.fasta_name}' \
      --mane \
      --canonical \
      --symbol \
      --biotype \
      --hgvs \
      --hgvsg \
      --shift_hgvs 1 \
      --numbers \
      --domains \
      --protein \
      --uniprot \
      --flag_pick \
      --pick_order mane_select,mane_plus_clinical,canonical,rank \
      --fork ${task.cpus} \
      --buffer_size ${params.vep_buffer_size} \
      --no_stats \
      --force_overwrite \
      "\${spliceai_args[@]}"

    printf 'vep_cache_dir\\t%s\\nvep_cache_version\\t%s\\nvep_buffer_size\\t%s\\nvep_forks\\t%s\\nspliceai\\t%s\\n' \
      "\${vep_cache_dir}" '${params.vep_cache_version}' '${params.vep_buffer_size}' \
      '${task.cpus}' "\${spliceai_status}" > vep-run-metadata.tsv
    """
}


process BUILD_REPORT {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    path annotated_vcf
    val sample_id

    output:
    path "${sample_id}.annotated.vcf.gz.tbi", emit: index
    path "${sample_id}.report.tsv", emit: report
    path 'report-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    bcftools index -t --threads ${task.cpus} '${annotated_vcf}'

    vep_variants=\$(bcftools index -n '${annotated_vcf}')
    vep_csq=\$(bcftools query -f '%INFO/CSQ\\n' '${annotated_vcf}' \
      | awk '\$0 != "." { n++ } END { print n+0 }')
    if (( vep_variants > 0 && vep_csq == 0 )); then
      echo 'ERROR: VEP produced no CSQ annotations' >&2
      exit 1
    fi

    printf 'CHROM\\tPOS\\tREF\\tALT\\tSYMBOL\\tTRANSCRIPT\\tCONSEQUENCE\\tIMPACT\\tHGVSc\\tHGVSp\\tMANE\\tPICK\\tGNOMAD_AF\\tGNOMAD_AF_MAX\\tCLNSIG\\tCLNSIGCONF\\tCLNDN\\tGT\\tDP\\tGQ\\tREASON_REPORTED\\n' \
      > '${sample_id}.report.tsv'

    bcftools +split-vep \
      -d \
      -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Feature\\t%Consequence\\t%IMPACT\\t%HGVSc\\t%HGVSp\\t%MANE_SELECT\\t%PICK\\t%INFO/gnomad_af\\t%INFO/gnomad_af_max\\t%INFO/CLNSIG\\t%INFO/CLNSIGCONF\\t%INFO/CLNDN\\t[%GT\\t%DP\\t%GQ]\\n' \
      -i 'IMPACT="HIGH" || IMPACT="MODERATE" || INFO/CLNSIG ~ "(^|[,/|])(Pathogenic|Likely_pathogenic)([,/|]|\$)" || INFO/CLNSIG ~ "Conflicting_classifications_of_pathogenicity" || INFO/CLNSIGCONF != "."' \
      '${annotated_vcf}' \
    | awk -F '\\t' -v OFS='\\t' -v metrics='report-metrics.tsv' \
        -v vep_variants="\${vep_variants}" -v vep_csq="\${vep_csq}" '
        {
          impact=(\$8 == "HIGH" || \$8 == "MODERATE")
          value=tolower(\$15)
          plp=(value ~ /(^|[,\\/|])(pathogenic|likely_pathogenic)([,\\/|]|\$)/)
          conflict=(value ~ /conflicting_classifications_of_pathogenicity/ || \$16 != ".")

          reason=""
          if (impact) reason="IMPACT"
          if (plp) reason=(reason == "" ? "P_LP" : reason "+P_LP")
          if (conflict) reason=(reason == "" ? "CONFLICT" : reason "+CONFLICT")

          print \$0, reason
          rows++
          if (plp) plp_rows++
          if (conflict) conflict_rows++
        }
        END {
          print "vep_variants", vep_variants > metrics
          print "vep_csq_annotated", vep_csq > metrics
          print "report_rows", rows+0 > metrics
          print "report_plp_rows", plp_rows+0 > metrics
          print "report_conflicting_rows", conflict_rows+0 > metrics
        }
      ' >> '${sample_id}.report.tsv'
    """
}


process BUILD_PROVENANCE {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 1
    memory '1 GB'
    time '30m'

    input:
    path counts
    path contig_filter_stats
    path report_metrics
    path vep_metadata
    val sample_id
    val reference_release
    val af_cutoff
    val af_max_cutoff

    output:
    path "${sample_id}.provenance.txt", emit: provenance

    script:
    """
    set -euo pipefail

    metric() {
      awk -F '\\t' -v key="\$1" '\$1 == key { print \$2; found=1; exit } END { if (!found) exit 1 }' "\$2"
    }

    variants_in=\$(metric variants_in '${counts}')
    variants_tiered=\$(metric variants_tiered '${counts}')
    variants_total=\$(metric variants_total '${contig_filter_stats}')
    variants_retained=\$(metric variants_retained '${contig_filter_stats}')
    variants_removed=\$(metric variants_removed_non_reference_contigs '${contig_filter_stats}')
    gnomad_af_annotated=\$(metric gnomad_af_annotated '${counts}')
    gnomad_af_max_annotated=\$(metric gnomad_af_max_annotated '${counts}')
    clinvar_annotated=\$(metric clinvar_annotated '${counts}')
    clinvar_plp=\$(metric clinvar_plp '${counts}')
    clinvar_conflicting=\$(metric clinvar_conflicting '${counts}')
    prefilter_plp_retained=\$(metric prefilter_plp_retained '${counts}')
    prefilter_conflicting_retained=\$(metric prefilter_conflicting_retained '${counts}')
    vep_csq_annotated=\$(metric vep_csq_annotated '${report_metrics}')
    report_rows=\$(metric report_rows '${report_metrics}')
    report_plp_rows=\$(metric report_plp_rows '${report_metrics}')
    report_conflicting_rows=\$(metric report_conflicting_rows '${report_metrics}')
    vep_cache_dir=\$(metric vep_cache_dir '${vep_metadata}')
    spliceai=\$(metric spliceai '${vep_metadata}')

    cat > '${sample_id}.provenance.txt' <<EOF
sample            ${sample_id}
date              \$(date -Iseconds)
assay_type        ${params.assay_type}
af_cutoff         ${af_cutoff}
af_max_cutoff     ${af_max_cutoff}
reference_release ${reference_release}
reference_mount   ${params.reference_dir}
reference         ${params.fasta_name}
clinvar           ${params.clinvar_name}
gnomad            ${params.gnomad_name}
vep               release ${params.vep_cache_version}, merged cache at \${vep_cache_dir}
mane              supplied by VEP release ${params.vep_cache_version} cache
spliceai          \${spliceai}
configured_images ${params.bcftools_image}
                  ${params.echtvar_image}
                  ${params.vep_image}
variants_in       \${variants_in}
variants_tiered   \${variants_tiered}
input_variants    \${variants_total}
primary_variants  \${variants_retained}
alt_removed       \${variants_removed}
gnomad_af_n       \${gnomad_af_annotated}
gnomad_af_max_n   \${gnomad_af_max_annotated}
clinvar_n         \${clinvar_annotated}
clinvar_plp_n     \${clinvar_plp}
clinvar_conflict_n \${clinvar_conflicting}
prefilter_plp_n   \${prefilter_plp_retained}
prefilter_conflict_n \${prefilter_conflicting_retained}
vep_csq_n         \${vep_csq_annotated}
report_rows       \${report_rows}
report_plp_rows   \${report_plp_rows}
report_conflict_rows \${report_conflicting_rows}

NOT VALIDATED FOR CLINICAL USE.
Missing: nhomalt (no recessive homozygote filter), internal AF panel,
SV/CNV/repeat annotation, constraint metrics, GIAB concordance run.
EOF
    """
}


process COLLECT_RESULTS {
    tag "${sample_id}"

    container params.bcftools_image
    cpus 1
    memory '1 GB'
    time '30m'

    publishDir params.dap_output_uri ?: params.output_dir, mode: 'copy', overwrite: false

    input:
    path annotated_vcf
    path annotated_index
    path report
    path provenance
    path counts
    path contig_filter_stats
    path report_metrics
    path vep_metadata
    val sample_id

    output:
    path 'results/*', emit: files

    script:
    """
    set -euo pipefail

    mkdir -p results
    cp '${annotated_vcf}' "results/${sample_id}.annotated.vcf.gz"
    cp '${annotated_index}' "results/${sample_id}.annotated.vcf.gz.tbi"
    cp '${report}' "results/${sample_id}.report.tsv"
    cp '${provenance}' "results/${sample_id}.provenance.txt"
    cp '${counts}' "results/${sample_id}.variant-counts.tsv"
    cp '${contig_filter_stats}' "results/${sample_id}.contig-filter.tsv"
    cp '${report_metrics}' "results/${sample_id}.report-metrics.tsv"
    cp '${vep_metadata}' "results/${sample_id}.vep-metadata.tsv"
    """
}


workflow {
    if (!params.vcf_uri || !params.vcf_index_uri) {
        error 'vcf_uri and vcf_index_uri are required'
    }
    if (!params.sample_id) {
        error 'sample_id is required'
    }
    if (!(params.sample_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) {
        error 'sample_id may contain only letters, numbers, period, underscore, and hyphen'
    }
    if (!params.reference_release) {
        error 'reference_release is required to make reference use cache-safe'
    }
    if (!(params.reference_release.toString() ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) {
        error 'reference_release may contain only letters, numbers, period, underscore, and hyphen'
    }
    if (!params.reference_dir || !params.reference_dir.toString().startsWith('/')) {
        error 'reference_dir must be an absolute mounted path'
    }
    if (!params.output_dir) {
        error 'output_dir is required'
    }
    if (params.vcf_sha256 && !(params.vcf_sha256 ==~ /[0-9a-fA-F]{64}/)) {
        error 'vcf_sha256 must contain 64 hexadecimal characters'
    }
    if (params.vcf_index_sha256 && !(params.vcf_index_sha256 ==~ /[0-9a-fA-F]{64}/)) {
        error 'vcf_index_sha256 must contain 64 hexadecimal characters'
    }
    if (!(params.assay_type in ['wgs', 'panel'])) {
        error 'assay_type must be wgs or panel'
    }
    if (params.vep_cache_subdir != 'auto' &&
        (!(params.vep_cache_subdir.toString() ==~ /[A-Za-z0-9][A-Za-z0-9._\\/-]*/) ||
         params.vep_cache_subdir.toString().contains('..') ||
         params.vep_cache_subdir.toString().startsWith('/'))) {
        error 'vep_cache_subdir must be auto or a safe relative path'
    }
    if (params.spliceai_indel_name &&
        (!(params.spliceai_indel_name.toString() ==~ /[A-Za-z0-9][A-Za-z0-9._\\/-]*/) ||
         params.spliceai_indel_name.toString().contains('..') ||
         params.spliceai_indel_name.toString().startsWith('/'))) {
        error 'spliceai_indel_name must be a safe relative path'
    }

    af_cutoff = params.af_cutoff as BigDecimal
    if (af_cutoff < 0 || af_cutoff > 1) {
        error 'af_cutoff must be between 0 and 1'
    }
    af_max_cutoff = params.af_max_cutoff as BigDecimal
    if (af_max_cutoff < 0 || af_max_cutoff > 1) {
        error 'af_max_cutoff must be between 0 and 1'
    }
    vep_buffer_size = params.vep_buffer_size as Integer
    if (vep_buffer_size < 1) {
        error 'vep_buffer_size must be a positive integer'
    }

    input_vcf = channel.fromPath(params.vcf_uri, checkIfExists: true)
    input_vcf_index = channel.fromPath(params.vcf_index_uri, checkIfExists: true)

    sample_id = channel.value(params.sample_id)
    vcf_sha256 = channel.value(params.vcf_sha256.toString().toLowerCase())
    vcf_index_sha256 = channel.value(params.vcf_index_sha256.toString().toLowerCase())
    reference_release = channel.value(params.reference_release)
    af_cutoff_value = channel.value(af_cutoff.toString())
    af_max_cutoff_value = channel.value(af_max_cutoff.toString())

    PRECHECK_REFERENCE(reference_release)

    FILTER_REFERENCE_CONTIGS(
        input_vcf,
        input_vcf_index,
        PRECHECK_REFERENCE.out.ready,
        sample_id,
        vcf_sha256,
        vcf_index_sha256,
        reference_release,
    )

    NORMALIZE(
        FILTER_REFERENCE_CONTIGS.out.vcf,
        FILTER_REFERENCE_CONTIGS.out.index,
        PRECHECK_REFERENCE.out.ready,
        sample_id,
        reference_release,
    )

    ANNOTATE_GNOMAD(
        NORMALIZE.out.vcf,
        NORMALIZE.out.index,
        PRECHECK_REFERENCE.out.ready,
        sample_id,
        reference_release,
    )

    ANNOTATE_CLINVAR(
        ANNOTATE_GNOMAD.out.vcf,
        PRECHECK_REFERENCE.out.ready,
        sample_id,
        reference_release,
    )

    PREFILTER(
        ANNOTATE_CLINVAR.out.vcf,
        ANNOTATE_CLINVAR.out.index,
        ANNOTATE_CLINVAR.out.metrics,
        sample_id,
        af_cutoff_value,
        af_max_cutoff_value,
    )

    ANNOTATE_VEP(
        PREFILTER.out.vcf,
        PREFILTER.out.index,
        PRECHECK_REFERENCE.out.ready,
        sample_id,
        reference_release,
    )

    BUILD_REPORT(
        ANNOTATE_VEP.out.vcf,
        sample_id,
    )

    BUILD_PROVENANCE(
        PREFILTER.out.counts,
        FILTER_REFERENCE_CONTIGS.out.stats,
        BUILD_REPORT.out.metrics,
        ANNOTATE_VEP.out.metadata,
        sample_id,
        reference_release,
        af_cutoff_value,
        af_max_cutoff_value,
    )

    COLLECT_RESULTS(
        ANNOTATE_VEP.out.vcf,
        BUILD_REPORT.out.index,
        BUILD_REPORT.out.report,
        BUILD_PROVENANCE.out.provenance,
        PREFILTER.out.counts,
        FILTER_REFERENCE_CONTIGS.out.stats,
        BUILD_REPORT.out.metrics,
        ANNOTATE_VEP.out.metadata,
        sample_id,
    )
}
