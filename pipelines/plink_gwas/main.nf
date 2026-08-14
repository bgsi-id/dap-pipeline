nextflow.enable.dsl = 2

include { READ_BUNDLES; WRITE_RESULTS } from '../base_nf/modules/bundle'

/*
 * Governed case/control GWAS, based on the Cloudfield GWASTutorial workflow:
 * source pVCF -> PLINK triples -> QC/PCA/KING -> Firth logistic GWAS -> plots.
 * The only phenotype authority is the two resolved immutable cohort packs.
 */
params.dap_input_manifest = null
params.dap_output_uri = null
params.cohort_id = null
params.output_dir = 'results'
params.n_pcs = 10
params.min_maf = 0.01
params.pre_qc_geno = 0.01
params.min_geno = 0.02
params.min_mind = 0.02
params.min_hwe = '1e-6'
params.het_sd = 0.1
params.king_cutoff = 0.0884
params.cpus = 8

process PREPARE_GWAS_INPUTS {
    tag "${params.cohort_id}"
    cpus 1
    memory '2 GB'
    container params.python_image

    input:
    val bundle

    output:
    path 'prepared'

    script:
    """
    set -euo pipefail
    mkdir prepared
    python3 - prepared '${groovy.json.JsonOutput.toJson(bundle).bytes.encodeBase64().toString()}' <<'PY'
import base64, json, pathlib, sys
bundle = json.loads(base64.b64decode(sys.argv[2])); out = pathlib.Path(sys.argv[1])
cohorts = bundle.get('cohorts', {})
control = cohorts.get('control', {}).get('samples', [])
patient = cohorts.get('patient', {}).get('samples', [])
if not control or not patient:
    raise SystemExit('control and patient cohorts must be non-empty')
if set(control) & set(patient):
    raise SystemExit('control and patient cohort packs overlap')
# PLINK case/control convention is 1=control and 2=case; zero is missing.
(out / 'keep.txt').write_text(''.join(f'{sample} {sample}\\n' for sample in control + patient))
(out / 'phenotype.txt').write_text('FID IID B1\\n' + ''.join(f'{sample} {sample} 1\\n' for sample in control) + ''.join(f'{sample} {sample} 2\\n' for sample in patient))
(out / 'provenance.json').write_text(json.dumps({
    'analysis_id': '${params.cohort_id}', 'manifest_sha256': bundle.get('manifest_sha256'),
    'dataset': bundle.get('dataset'), 'control': cohorts['control'], 'patient': cohorts['patient'],
}, indent=2, sort_keys=True) + '\\n')
PY
    """
}

process CONVERT_TO_PLINK {
    tag "${params.cohort_id}"
    cpus params.cpus
    memory '16 GB'
    container params.plink_image

    input:
    path prepared
    tuple path(pvcf), path(pvcf_index)

    output:
    path 'raw.*', emit: raw
    path 'prepared', emit: prepared

    script:
    """
    set -euo pipefail
    plink --vcf '${pvcf}' --id-delim '_' --keep prepared/keep.txt --keep-allele-order --make-bed --out raw
    awk '{ print \$2 }' raw.fam | sort > imported-samples.txt
    awk '{ print \$2 }' prepared/keep.txt | sort > requested-samples.txt
    if ! diff -u requested-samples.txt imported-samples.txt; then
      echo 'Resolved cohort membership is not fully present in the pVCF' >&2
      exit 1
    fi
    """
}

process PRE_GWAS_QC {
    tag "${params.cohort_id}"
    cpus params.cpus
    memory '24 GB'
    container params.plink_image

    input:
    path raw
    path prepared

    output:
    path 'pregwas'

    script:
    """
    set -euo pipefail
    mkdir pregwas
    # Cloudfield pre-GWAS QC: missingness, MAF/HWE, heterozygosity, LD, KING, PCA.
    plink --bfile raw --missing --freq --hardy --out pregwas/basic
    plink --bfile raw --maf ${params.min_maf} --geno ${params.pre_qc_geno} --mind ${params.min_mind} --hwe ${params.min_hwe} --indep-pairwise 50 5 0.2 --out pregwas/qc-prune
    plink --bfile raw --extract pregwas/qc-prune.prune.in --het --out pregwas/qc-het
    awk 'NR > 1 && (\$6 > ${params.het_sd} || \$6 < -${params.het_sd}) { print \$1, \$2 }' pregwas/qc-het.het > pregwas/high-het.sample
    plink --bfile raw --geno ${params.min_geno} --mind ${params.min_mind} --hwe ${params.min_hwe} --remove pregwas/high-het.sample --keep-allele-order --make-bed --out pregwas/clean
    """
}

process PCA_AND_KING {
    tag "${params.cohort_id}"
    cpus params.cpus
    memory '24 GB'
    container params.plink2_image

    input:
    path pregwas

    output:
    path 'pca'

    script:
    """
    set -euo pipefail
    mkdir pca
    plink2 --bfile pregwas/clean --maf ${params.min_maf} --indep-pairwise 500 50 0.2 --out pca/pca-prune
    plink2 --bfile pregwas/clean --extract pca/pca-prune.prune.in --king-cutoff ${params.king_cutoff} --out pca/king
    plink2 --bfile pregwas/clean --keep pca/king.king.cutoff.in.id --extract pca/pca-prune.prune.in --freq counts --pca allele-wts ${params.n_pcs} --out pca/pca
    plink2 --bfile pregwas/clean --read-freq pca/pca.acount --score pca/pca.eigenvec.allele 2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-${params.n_pcs + 6} --out pca/pca-projected
    """
}

process RUN_GWAS {
    tag "${params.cohort_id}"
    cpus params.cpus
    memory '32 GB'
    container params.plink2_image

    input:
    path pregwas
    path pca
    path prepared

    output:
    path 'association'

    script:
    def pcNames = (1..params.n_pcs).collect { "PC${it}_AVG" }.join(',')
    """
    set -euo pipefail
    mkdir association
    plink2 --bfile pregwas/clean --keep pca/king.king.cutoff.in.id --pheno prepared/phenotype.txt --pheno-name B1 --maf ${params.min_maf} --covar pca/pca-projected.sscore --covar-name ${pcNames} --glm hide-covar firth firth-residualize single-prec-cc --threads ${task.cpus} --out association/${params.cohort_id}
    """
}

process POST_GWAS {
    tag "${params.cohort_id}"
    cpus 1
    memory '4 GB'
    container params.plot_image

    input:
    path association
    path pca
    path pregwas
    path prepared

    output:
    path "${params.output_dir}"

    script:
    """
    set -euo pipefail
    mkdir -p '${params.output_dir}'
    cp association/${params.cohort_id}.B1.glm.firth '${params.output_dir}/${params.cohort_id}.sumstats.tsv'
    cp prepared/provenance.json '${params.output_dir}/${params.cohort_id}.provenance.json'
    cp pregwas/basic.{imiss,lmiss,frq,hwe} '${params.output_dir}/'
    cp pregwas/{qc-prune.prune.in,high-het.sample} '${params.output_dir}/'
    cp pca/{king.king.cutoff.in.id,pca.eigenvec,pca-projected.sscore} '${params.output_dir}/'
    python3 - association/${params.cohort_id}.B1.glm.firth '${params.output_dir}' '${params.cohort_id}' <<'PY'
import csv, json, math, pathlib, sys
source, out, prefix = map(pathlib.Path, sys.argv[1:])
rows = []
with source.open() as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        try: p = float(row['P'])
        except (ValueError, KeyError): continue
        if 0 < p <= 1: rows.append((row, p))
p_ranked = sorted(rows, key=lambda item: item[1])
with (out / f'{prefix}.top-hits.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(p_ranked[0][0]) if p_ranked else ['ID','P'], delimiter='\t')
    writer.writeheader(); writer.writerows(row for row, _ in p_ranked[:100])
summary = {'tested_variants': len(rows), 'top_hit': p_ranked[0][0] if p_ranked else None}
(out / f'{prefix}.summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True))
try:
    import matplotlib.pyplot as plt

    def chromosome_label(value):
        label = str(value)
        return label[3:] if label.lower().startswith('chr') else label

    def chromosome_key(value):
        label = chromosome_label(value).upper()
        aliases = {'X': 23, 'Y': 24, 'XY': 25, 'M': 26, 'MT': 26}
        try:
            number = int(label)
            return (number if number > 0 else 10_000, label)
        except ValueError:
            return (aliases.get(label, 10_000), label)

    genomic_rows = []
    for row, p in rows:
        chromosome = str(row.get('#CHROM') or row.get('CHROM') or '').strip()
        try:
            position = int(row.get('POS', ''))
        except (TypeError, ValueError):
            continue
        if chromosome and position >= 0:
            genomic_rows.append((chromosome_key(chromosome), chromosome, position, p))
    genomic_rows.sort(key=lambda item: (item[0], item[2]))

    if genomic_rows:
        chromosome_max = {}
        chromosome_labels = {}
        for key, label, position, _ in genomic_rows:
            chromosome_max[key] = max(chromosome_max.get(key, 0), position)
            chromosome_labels.setdefault(key, chromosome_label(label))

        offsets = {}
        cursor = 0
        for key in sorted(chromosome_max):
            offsets[key] = cursor
            cursor += chromosome_max[key] + 1

        x = [offsets[key] + position for key, _, position, _ in genomic_rows]
        y = [-math.log10(p) for _, _, _, p in genomic_rows]
        chromosome_index = {key: index for index, key in enumerate(sorted(chromosome_max))}
        colors = ['#3155b7' if chromosome_index[key] % 2 == 0 else '#b21f35' for key, _, _, _ in genomic_rows]
        ticks = [offsets[key] + chromosome_max[key] / 2 for key in sorted(chromosome_max)]
        labels = [chromosome_labels[key] for key in sorted(chromosome_max)]

        plt.figure(figsize=(16, 7))
        plt.scatter(x, y, c=colors, s=3, linewidths=0, rasterized=True)
        plt.axhline(-math.log10(5e-8), color='#b21f35', linestyle='--', linewidth=1, label='P = 5e-8')
        plt.axhline(-math.log10(1e-5), color='#777777', linestyle=':', linewidth=1, label='P = 1e-5')
        plt.xticks(ticks, labels)
        plt.xlabel('Chromosome')
        plt.ylabel('-log10(P)')
        plt.title('Manhattan plot')
        plt.legend(frameon=False, loc='upper right')
        plt.tight_layout()
        plt.savefig(out / f'{prefix}.manhattan.png', dpi=180)
        plt.close()

    if p_ranked:
        observed = [-math.log10(p) for _, p in p_ranked]
        expected = [-math.log10((i + 0.5) / len(p_ranked)) for i in range(len(p_ranked))]
        limit = max(expected[0], observed[0])
        plt.figure(figsize=(7, 7))
        plt.scatter(expected, observed, s=3, linewidths=0, rasterized=True)
        plt.plot([0, limit], [0, limit], color='#b21f35', linewidth=1)
        plt.xlim(0, limit * 1.02)
        plt.ylim(0, limit * 1.02)
        plt.xlabel('Expected -log10(P)')
        plt.ylabel('Observed -log10(P)')
        plt.title('QQ plot')
        plt.tight_layout()
        plt.savefig(out / f'{prefix}.qq.png', dpi=180)
        plt.close()
except ImportError:
    pass
PY
    """
}

workflow {
    if (!params.dap_input_manifest || !params.cohort_id) error 'dap_input_manifest and cohort_id are required'
    if (!params.dap_output_uri) error 'dap_output_uri is required for project output publishing'
    bundles = READ_BUNDLES(Channel.value(file(params.dap_input_manifest)))
    release_inputs = bundles.bundle.map { bundle ->
        def pvcf = bundle.release_assets.pvcf
        def pvcf_index = bundle.release_assets.pvcf_index
        if (!pvcf || !pvcf_index) error 'PLINK GWAS requires pvcf and pvcf_index release assets'
        tuple(file(pvcf.access_uri), file(pvcf_index.access_uri))
    }
    prepared = PREPARE_GWAS_INPUTS(bundles.bundle)
    converted = CONVERT_TO_PLINK(prepared, release_inputs)
    qc = PRE_GWAS_QC(converted.raw, converted.prepared)
    pca = PCA_AND_KING(qc)
    association = RUN_GWAS(qc, pca, converted.prepared)
    results = POST_GWAS(association, pca, qc, converted.prepared)
    WRITE_RESULTS(results)
}
