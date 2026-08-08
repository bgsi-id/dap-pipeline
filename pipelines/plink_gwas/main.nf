nextflow.enable.dsl = 2

/*
 * Governed case/control GWAS, based on the Cloudfield GWASTutorial workflow:
 * source pVCF -> PLINK triples -> QC/PCA/KING -> Firth logistic GWAS -> plots.
 * The only phenotype authority is the two resolved immutable cohort packs.
 */
params.dap_input_manifest = null
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
    path manifest

    output:
    path 'prepared'

    script:
    """
    set -euo pipefail
    mkdir prepared
    python3 - '${manifest}' prepared <<'PY'
import json, pathlib, sys
manifest = json.load(open(sys.argv[1])); out = pathlib.Path(sys.argv[2])
if manifest.get('schema') != 'urn:bgsi:dap:resolved-inputs:2':
    raise SystemExit('requires dap resolved-inputs v2')
cohorts = manifest.get('cohorts', {})
control = cohorts.get('control', {}).get('samples', [])
patient = cohorts.get('patient', {}).get('samples', [])
if not control or not patient:
    raise SystemExit('control and patient cohorts must be non-empty')
if set(control) & set(patient):
    raise SystemExit('control and patient cohort packs overlap')
assets = {item['role']: item for item in manifest.get('release_assets', [])}
for role in ('pvcf', 'pvcf_index'):
    if not assets.get(role, {}).get('access_uri'):
        raise SystemExit(f'missing release asset: {role}')
# PLINK case/control convention is 1=control and 2=case; zero is missing.
(out / 'keep.txt').write_text(''.join(f'{sample} {sample}\\n' for sample in control + patient))
(out / 'phenotype.txt').write_text('FID IID B1\\n' + ''.join(f'{sample} {sample} 1\\n' for sample in control) + ''.join(f'{sample} {sample} 2\\n' for sample in patient))
(out / 'provenance.json').write_text(json.dumps({
    'analysis_id': '${params.cohort_id}', 'manifest_sha256': manifest.get('manifest_sha256'),
    'dataset': manifest.get('dataset'), 'control': cohorts['control'], 'patient': cohorts['patient'],
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
    path pvcf
    path pvcf_index

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
    plink2 --bfile pregwas/clean --keep pca/king.king.cutoff.in.id --extract pca/pca-prune.prune.in --freq counts --pca approx allele-wts ${params.n_pcs} --out pca/pca
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
rows.sort(key=lambda item: item[1])
with (out / f'{prefix}.top-hits.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0][0]) if rows else ['ID','P'], delimiter='\t')
    writer.writeheader(); writer.writerows(row for row, _ in rows[:100])
summary = {'tested_variants': len(rows), 'top_hit': rows[0][0] if rows else None}
(out / f'{prefix}.summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
try:
    import matplotlib.pyplot as plt
    chrom = [str(row.get('#CHROM') or row.get('CHROM') or '') for row, _ in rows]
    y = [-math.log10(p) for _, p in rows]
    x = list(range(len(rows)))
    plt.figure(figsize=(16, 7)); plt.scatter(x, y, c=[hash(c) % 2 for c in chrom], s=2, cmap='coolwarm'); plt.xlabel('Variant'); plt.ylabel('-log10(P)'); plt.title('Manhattan plot'); plt.tight_layout(); plt.savefig(out / f'{prefix}.manhattan.png', dpi=180); plt.close()
    expected = [-math.log10((i + .5) / len(rows)) for i in range(len(rows))]
    observed = sorted(y)
    plt.figure(figsize=(7, 7)); plt.scatter(expected, observed, s=3); plt.plot([0, max(expected, default=1)], [0, max(expected, default=1)], 'r-'); plt.xlabel('Expected -log10(P)'); plt.ylabel('Observed -log10(P)'); plt.title('QQ plot'); plt.tight_layout(); plt.savefig(out / f'{prefix}.qq.png', dpi=180); plt.close()
except ImportError:
    pass
PY
    """
}

workflow {
    if (!params.dap_input_manifest || !params.cohort_id) error 'dap_input_manifest and cohort_id are required'
    def resolved = new groovy.json.JsonSlurper().parse(file(params.dap_input_manifest))
    def assets = (resolved.release_assets ?: []).collectEntries { asset ->
        [(asset.role.toString().toLowerCase()): asset]
    }
    def pvcf_uri = assets.pvcf?.access_uri
    def pvcf_index_uri = assets.pvcf_index?.access_uri
    if (!pvcf_uri || !pvcf_index_uri) error 'resolved inputs must include pvcf and pvcf_index release assets'
    prepared = PREPARE_GWAS_INPUTS(file(params.dap_input_manifest))
    converted = CONVERT_TO_PLINK(prepared, file(pvcf_uri), file(pvcf_index_uri))
    qc = PRE_GWAS_QC(converted.raw, converted.prepared)
    pca = PCA_AND_KING(qc)
    association = RUN_GWAS(qc, pca, converted.prepared)
    POST_GWAS(association, pca, qc, converted.prepared)
}
