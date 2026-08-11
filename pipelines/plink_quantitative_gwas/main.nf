nextflow.enable.dsl = 2

include { WRITE_RESULTS } from '../base_nf/modules/bundle'

params.genotype_uri = null
params.genotype_format = 'vcf'
params.bgen_sample_uri = ''
params.phenotype_uri = null
params.dap_output_uri = null
params.cohort_id = null
params.sample_id_column = 'IID'
params.phenotype_column = null
params.output_dir = 'results'
params.n_pcs = 10
params.min_maf = 0.01
params.min_mac = 20
params.min_geno = 0.02
params.min_mind = 0.02
params.min_hwe = '1e-6'
params.king_cutoff = 0.0884
params.cpus = 8
params.plink_memory = '32 GB'

process PREPARE_PHENOTYPE {
    tag params.cohort_id
    cpus 1
    memory '2 GB'
    container params.python_image

    input:
    path phenotypes

    output:
    path 'prepared'

    script:
    """
    set -euo pipefail
    mkdir prepared
    python3 - '${phenotypes}' prepared/phenotype.tsv prepared/phenotype-metrics.tsv <<'PY'
import csv, math, pathlib, statistics, sys
source, output, metrics = map(pathlib.Path, sys.argv[1:])
with source.open(newline='') as handle:
    sample = handle.read(4096); handle.seek(0)
    dialect = csv.Sniffer().sniff(sample, delimiters='\t,')
    rows = list(csv.DictReader(handle, dialect=dialect))
sid, trait = '${params.sample_id_column}', '${params.phenotype_column}'
if not rows or sid not in rows[0] or trait not in rows[0]:
    raise SystemExit(f'phenotype file must contain {sid!r} and {trait!r}')
seen, values, valid = set(), [], []
for row in rows:
    sample_id = (row.get(sid) or '').strip()
    if not sample_id or sample_id in seen:
        raise SystemExit('sample IDs must be non-empty and unique')
    seen.add(sample_id)
    raw = (row.get(trait) or '').strip()
    if raw.lower() in {'', 'na', 'nan', '.'}: continue
    try: value = float(raw)
    except ValueError: raise SystemExit(f'quantitative phenotype is not numeric for {sample_id}')
    if not math.isfinite(value): continue
    values.append(value); valid.append((sample_id, value))
if len(valid) < 20 or len(set(values)) < 3:
    raise SystemExit('quantitative GWAS requires at least 20 samples and 3 distinct values')
with output.open('w') as handle:
    handle.write('#FID\tIID\tPHENO1\n')
    for sample_id, value in valid: handle.write(f'{sample_id}\t{sample_id}\t{value:.12g}\n')
with metrics.open('w') as handle:
    handle.write(f'samples_in\t{len(rows)}\n')
    handle.write(f'samples_with_phenotype\t{len(valid)}\n')
    handle.write(f'phenotype_mean\t{statistics.fmean(values):.8g}\n')
    handle.write(f'phenotype_sd\t{statistics.stdev(values):.8g}\n')
PY
    """

    stub:
    """
    mkdir prepared
    printf '#FID\tIID\tPHENO1\nS1\tS1\t1.0\n' > prepared/phenotype.tsv
    printf 'samples_in\t20\nsamples_with_phenotype\t20\n' > prepared/phenotype-metrics.tsv
    """
}

process IMPORT_AND_QC {
    tag params.cohort_id
    cpus params.cpus
    memory params.plink_memory
    container params.plink2_image

    input:
    path genotype
    path bgen_sample
    path prepared

    output:
    path 'qc', emit: qc
    path 'prepared', emit: prepared

    script:
    def importArgs = params.genotype_format == 'bgen'
        ? "--bgen '${genotype}' ref-first --sample '${bgen_sample}'"
        : "--vcf '${genotype}'"
    """
    set -euo pipefail
    mkdir qc
    plink2 ${importArgs} --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 missing --make-pgen --out qc/raw
    plink2 --pfile qc/raw --pheno prepared/phenotype.tsv --pheno-name PHENO1 \
      --geno ${params.min_geno} --mind ${params.min_mind} --maf ${params.min_maf} \
      --hwe ${params.min_hwe} midp keep-fewhet --make-pgen --out qc/clean
    plink2 --pfile qc/clean --missing sample-only --freq --hardy --out qc/metrics
    """

    stub:
    """
    mkdir qc
    touch qc/clean.pgen qc/clean.pvar qc/clean.psam qc/metrics.smiss qc/metrics.afreq qc/metrics.hardy
    """
}

process PCA_AND_RELATEDNESS {
    tag params.cohort_id
    cpus params.cpus
    memory params.plink_memory
    container params.plink2_image

    input:
    path qc

    output:
    path 'structure'

    script:
    """
    set -euo pipefail
    mkdir structure
    plink2 --pfile qc/clean --indep-pairwise 500kb 0.2 --out structure/prune
    plink2 --pfile qc/clean --extract structure/prune.prune.in --king-cutoff ${params.king_cutoff} --out structure/king
    plink2 --pfile qc/clean --keep structure/king.king.cutoff.in.id --extract structure/prune.prune.in \
      --pca approx ${params.n_pcs} --out structure/pca
    """

    stub:
    """
    mkdir structure
    printf '#FID\tIID\nS1\tS1\n' > structure/king.king.cutoff.in.id
    printf '#FID\tIID\tPC1\nS1\tS1\t0.1\n' > structure/pca.eigenvec
    """
}

process RUN_LINEAR_GWAS {
    tag params.cohort_id
    cpus params.cpus
    memory params.plink_memory
    container params.plink2_image

    input:
    path qc
    path structure
    path prepared

    output:
    path 'association'

    script:
    def pcs = (1..params.n_pcs).collect { "PC${it}" }.join(',')
    """
    set -euo pipefail
    mkdir association
    plink2 --pfile qc/clean --keep structure/king.king.cutoff.in.id \
      --pheno prepared/phenotype.tsv --pheno-name PHENO1 \
      --covar structure/pca.eigenvec --covar-name ${pcs} \
      --mac ${params.min_mac} --glm hide-covar qt-residualize \
      --threads ${task.cpus} --out association/${params.cohort_id}
    test -s association/${params.cohort_id}.PHENO1.glm.linear
    """

    stub:
    """
    mkdir association
    printf '#CHROM\tPOS\tID\tREF\tALT\tP\n1\t1\tv1\tA\tG\t0.5\n' > association/${params.cohort_id}.PHENO1.glm.linear
    """
}

process COLLECT_RESULTS {
    tag params.cohort_id
    cpus 1
    memory '4 GB'
    container params.python_image

    input:
    path association
    path qc
    path structure
    path prepared

    output:
    path 'results'

    script:
    """
    set -euo pipefail
    mkdir results
    cp association/${params.cohort_id}.PHENO1.glm.linear results/${params.cohort_id}.sumstats.tsv
    cp structure/pca.eigenvec results/${params.cohort_id}.pcs.tsv
    cp prepared/phenotype-metrics.tsv results/${params.cohort_id}.phenotype-metrics.tsv
    cp qc/metrics.* results/
    python3 - results/${params.cohort_id}.sumstats.tsv results/${params.cohort_id}.top-hits.tsv <<'PY'
import csv, pathlib, sys
source, output = map(pathlib.Path, sys.argv[1:]); rows=[]
with source.open() as handle:
    reader=csv.DictReader(handle, delimiter='\t')
    fields=reader.fieldnames or []
    for row in reader:
        try: p=float(row['P'])
        except (KeyError, ValueError, TypeError): continue
        if 0 < p <= 1: rows.append((p,row))
rows.sort(key=lambda item:item[0])
with output.open('w') as handle:
    writer=csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader()
    writer.writerows(row for _,row in rows[:100])
PY
    printf 'pipeline\tGWAS (PLINK2) - quantitative traits\nphenotype\t%s\nmodel\tlinear --glm with PCA covariates\n' '${params.phenotype_column}' > results/${params.cohort_id}.provenance.tsv
    """

    stub:
    """
    mkdir results
    touch results/${params.cohort_id}.sumstats.tsv results/${params.cohort_id}.top-hits.tsv results/${params.cohort_id}.pcs.tsv results/${params.cohort_id}.provenance.tsv
    """
}

workflow {
    if (!params.genotype_uri || !params.phenotype_uri || !params.phenotype_column || !params.cohort_id) {
        error 'genotype_uri, phenotype_uri, phenotype_column, and cohort_id are required'
    }
    if (!(params.genotype_format in ['vcf', 'bgen'])) error 'genotype_format must be vcf or bgen'
    if (params.genotype_format == 'bgen' && !params.bgen_sample_uri) error 'bgen_sample_uri is required for BGEN input'
    if (!(params.cohort_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) error 'invalid cohort_id'
    column_pattern = /[A-Za-z_][A-Za-z0-9_.]{0,63}/
    if (!(params.sample_id_column ==~ column_pattern) || !(params.phenotype_column ==~ column_pattern)) error 'invalid phenotype column name'

    genotype = channel.fromPath(params.genotype_uri, checkIfExists: true)
    phenotypes = channel.fromPath(params.phenotype_uri, checkIfExists: true)
    bgen_sample = params.bgen_sample_uri ? channel.fromPath(params.bgen_sample_uri, checkIfExists: true) : channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)
    prepared = PREPARE_PHENOTYPE(phenotypes)
    imported = IMPORT_AND_QC(genotype, bgen_sample, prepared)
    structure = PCA_AND_RELATEDNESS(imported.qc)
    association = RUN_LINEAR_GWAS(imported.qc, structure, imported.prepared)
    results = COLLECT_RESULTS(association, imported.qc, structure, imported.prepared)
    WRITE_RESULTS(results)
}
