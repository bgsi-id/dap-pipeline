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
params.af_cutoff = 0.01

params.fasta_name = 'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
params.gnomad_name = 'gnomad_v4.1.zip'
params.clinvar_name = 'clinvar.chr.vcf.gz'
params.spliceai_name = 'spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz'
params.vep_cache_subdir = 'vep'
params.vep_cache_version = '116'

params.bcftools_image = 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
params.echtvar_image = 'quay.io/biocontainers/echtvar:0.2.2--h4349ce8_0'
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

    test -s '${params.reference_dir}/${params.fasta_name}'
    test -s '${params.reference_dir}/${params.fasta_name}.fai'
    test -s '${params.reference_dir}/${params.gnomad_name}'
    test -s '${params.reference_dir}/${params.clinvar_name}'
    test -s '${params.reference_dir}/${params.clinvar_name}.tbi'
    test -d '${params.reference_dir}/${params.vep_cache_subdir}'

    printf '%s\\n' '${reference_release}' > reference.ready
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

    script:
    """
    set -euo pipefail

    bcftools index -t --threads ${task.cpus} '${frequency_vcf}'
    bcftools annotate \
      -a '${params.reference_dir}/${params.clinvar_name}' \
      -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVI \
      --threads ${task.cpus} \
      -Oz \
      -o annot.vcf.gz \
      '${frequency_vcf}'
    bcftools index -t --threads ${task.cpus} annot.vcf.gz
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
    val sample_id
    val af_cutoff

    output:
    path 'tiered.vcf.gz', emit: vcf
    path 'tiered.vcf.gz.tbi', emit: index
    path 'variant-counts.tsv', emit: counts

    script:
    """
    set -euo pipefail

    bcftools filter \
      -i '((INFO/gnomad_af < ${af_cutoff} || INFO/gnomad_af = -1) && (INFO/gnomad_af_max < 0.02 || INFO/gnomad_af_max = -1)) || INFO/CLNSIG ~ "athogenic"' \
      --threads ${task.cpus} \
      -Oz \
      -o tiered.vcf.gz \
      '${annotated_vcf}'
    bcftools index -t --threads ${task.cpus} tiered.vcf.gz

    variants_in=\$(bcftools index -n '${annotated_vcf}')
    variants_tiered=\$(bcftools index -n tiered.vcf.gz)
    printf 'variants_in\\t%s\\nvariants_tiered\\t%s\\n' \
      "\${variants_in}" "\${variants_tiered}" > variant-counts.tsv
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

    script:
    """
    set -euo pipefail

    spliceai_args=()
    if [[ -s '${params.reference_dir}/${params.spliceai_name}' && -s '${params.reference_dir}/${params.spliceai_name}.tbi' ]]; then
      spliceai_args+=(--plugin 'SpliceAI,snv=${params.reference_dir}/${params.spliceai_name}')
    fi

    vep \
      --input_file '${tiered_vcf}' \
      --output_file '${sample_id}.annotated.vcf.gz' \
      --vcf \
      --compress_output bgzip \
      --offline \
      --cache \
      --dir_cache '${params.reference_dir}/${params.vep_cache_subdir}' \
      --merged \
      --cache_version '${params.vep_cache_version}' \
      --assembly GRCh38 \
      --fasta '${params.reference_dir}/${params.fasta_name}' \
      --mane \
      --mane_select \
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
      --pick_order mane_select,mane_plus_clinical,canonical,rank \
      --fork ${task.cpus} \
      --buffer_size 50000 \
      --no_stats \
      --force_overwrite \
      "\${spliceai_args[@]}"
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

    script:
    """
    set -euo pipefail

    bcftools index -t --threads ${task.cpus} '${annotated_vcf}'
    {
      printf 'CHROM\\tPOS\\tREF\\tALT\\tSYMBOL\\tTRANSCRIPT\\tCONSEQUENCE\\tIMPACT\\tHGVSc\\tHGVSp\\tMANE\\tGNOMAD_AF\\tGNOMAD_AF_MAX\\tCLNSIG\\tCLNDN\\tGT\\tDP\\tGQ\\n'
      bcftools +split-vep \
        -d \
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Feature\\t%Consequence\\t%IMPACT\\t%HGVSc\\t%HGVSp\\t%MANE_SELECT\\t%INFO/gnomad_af\\t%INFO/gnomad_af_max\\t%INFO/CLNSIG\\t%INFO/CLNDN\\t[%GT\\t%DP\\t%GQ]\\n' \
        -i 'IMPACT="HIGH" || IMPACT="MODERATE" || CLNSIG ~ "athogenic"' \
        '${annotated_vcf}'
    } > '${sample_id}.report.tsv'
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
    val sample_id
    val reference_release
    val af_cutoff

    output:
    path "${sample_id}.provenance.txt", emit: provenance

    script:
    """
    set -euo pipefail

    variants_in=\$(sed -n '1s/^[^[:space:]]*[[:space:]]*//p' '${counts}')
    variants_tiered=\$(sed -n '2s/^[^[:space:]]*[[:space:]]*//p' '${counts}')
    variants_total=\$(sed -n '1s/^[^[:space:]]*[[:space:]]*//p' '${contig_filter_stats}')
    variants_retained=\$(sed -n '2s/^[^[:space:]]*[[:space:]]*//p' '${contig_filter_stats}')
    variants_removed=\$(sed -n '3s/^[^[:space:]]*[[:space:]]*//p' '${contig_filter_stats}')
    spliceai='NOT APPLIED'
    if [[ -s '${params.reference_dir}/${params.spliceai_name}' && -s '${params.reference_dir}/${params.spliceai_name}.tbi' ]]; then
      spliceai='ensembl_mane_v1.4'
    fi

    cat > '${sample_id}.provenance.txt' <<EOF
sample            ${sample_id}
date              \$(date -Iseconds)
af_cutoff         ${af_cutoff}
reference_release ${reference_release}
reference_mount   ${params.reference_dir}
reference         ${params.fasta_name}
vep               release ${params.vep_cache_version}, merged cache
mane              v1.5 summary (report layer)
gnomad            v4.1 genomes, AF and AF_grpmax only
spliceai          \${spliceai}
images            ${params.bcftools_image}
                  ${params.echtvar_image}
                  ${params.vep_image}
variants_in       \${variants_in}
variants_tiered   \${variants_tiered}
input_variants    \${variants_total}
primary_variants  \${variants_retained}
alt_removed       \${variants_removed}

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

    publishDir params.output_dir, mode: 'copy', overwrite: false

    input:
    path annotated_vcf
    path annotated_index
    path report
    path provenance
    path counts
    path contig_filter_stats
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

    af_cutoff = params.af_cutoff as BigDecimal
    if (af_cutoff < 0 || af_cutoff > 1) {
        error 'af_cutoff must be between 0 and 1'
    }

    input_vcf = channel.fromPath(params.vcf_uri, checkIfExists: true)
    input_vcf_index = channel.fromPath(params.vcf_index_uri, checkIfExists: true)

    sample_id = channel.value(params.sample_id)
    vcf_sha256 = channel.value(params.vcf_sha256.toString().toLowerCase())
    vcf_index_sha256 = channel.value(params.vcf_index_sha256.toString().toLowerCase())
    reference_release = channel.value(params.reference_release)
    af_cutoff_value = channel.value(af_cutoff.toString())

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
        sample_id,
        af_cutoff_value,
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
        sample_id,
        reference_release,
        af_cutoff_value,
    )

    COLLECT_RESULTS(
        ANNOTATE_VEP.out.vcf,
        BUILD_REPORT.out.index,
        BUILD_REPORT.out.report,
        BUILD_PROVENANCE.out.provenance,
        PREFILTER.out.counts,
        FILTER_REFERENCE_CONTIGS.out.stats,
        sample_id,
    )
}
