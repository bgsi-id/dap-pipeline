nextflow.enable.dsl = 2

/*
 * Cohort GWAS (Hail)
 *
 * import -> genotype/variant/sample QC -> PCA -> association -> Manhattan/QQ
 *
 * Sized for 500-2000 sample cohorts running Hail in Spark local mode inside a
 * single container. Every stage emits keyed metrics and terminates rather than
 * publishing an association result derived from an unusable cohort.
 */

params.genotype_uri = null
params.genotype_format = 'mt'
params.genotype_sha256 = ''
params.bgen_sample_uri = ''
params.phenotype_uri = null
params.phenotype_sha256 = ''

params.cohort_id = null
params.output_dir = 'results'
params.reference_genome = 'GRCh38'

params.sample_id_column = 's'
params.phenotype_column = null
params.phenotype_type = 'quantitative'
params.covariate_columns = ''

params.min_gq = 20
params.min_dp = 10
params.min_call_rate_variant = 0.95
params.min_call_rate_sample = 0.95
params.min_maf = 0.01
params.min_hwe_p = 1e-6

params.n_pcs = 10
params.ld_prune = true
params.ld_prune_r2 = 0.2
params.ld_prune_window = 500000
params.pca_max_variants = 200000

params.min_samples_retained = 100
params.min_variants_retained = 10000
params.min_lambda_gc = 0.8
params.max_lambda_gc = 1.2

params.hail_image = 'hailgenetics/hail:0.2.133'
params.hail_driver_memory = '24g'


process IMPORT_GENOTYPES {
    tag "${cohort_id}"

    container params.hail_image
    cpus 8
    memory '32 GB'
    time '8h'

    input:
    path genotype
    path bgen_sample
    val cohort_id
    val genotype_sha256

    output:
    path 'cohort.mt', emit: mt
    path 'import-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    if [[ -n '${genotype_sha256}' && -f '${genotype}' ]]; then
      echo '${genotype_sha256}  ${genotype}' | sha256sum -c -
    fi

    mkdir -p hail-tmp

    python3 <<'PY'
import hail as hl

hl.init(
    master='local[${task.cpus}]',
    tmp_dir='hail-tmp',
    local_tmpdir='hail-tmp',
    default_reference='${params.reference_genome}',
    spark_conf={'spark.driver.memory': '${params.hail_driver_memory}'},
    log='hail-import.log',
)

fmt = '${params.genotype_format}'
if fmt == 'mt':
    mt = hl.read_matrix_table('${genotype}')
elif fmt == 'vcf':
    mt = hl.import_vcf(
        '${genotype}',
        force_bgz=True,
        reference_genome='${params.reference_genome}',
        array_elements_required=False,
    )
elif fmt == 'bgen':
    hl.index_bgen('${genotype}', reference_genome='${params.reference_genome}')
    sample_file = '${bgen_sample}'
    if sample_file.endswith('NO_FILE'):
        mt = hl.import_bgen('${genotype}', entry_fields=['GT', 'dosage'])
    else:
        mt = hl.import_bgen(
            '${genotype}',
            entry_fields=['GT', 'dosage'],
            sample_file=sample_file,
        )
else:
    raise SystemExit(f'unsupported genotype_format: {fmt}')

if 'GT' not in mt.entry:
    raise SystemExit('input has no GT entry field; GWAS requires hard calls')

mt = mt.key_rows_by('locus', 'alleles')

n_variants, n_samples = mt.count()
if n_samples == 0:
    raise SystemExit('input contains no samples')
if n_variants == 0:
    raise SystemExit('input contains no variants')

mt.write('cohort.mt', overwrite=True)

with open('import-metrics.tsv', 'w') as handle:
    handle.write(f'genotype_format\\t{fmt}\\n')
    handle.write(f'samples_imported\\t{n_samples}\\n')
    handle.write(f'variants_imported\\t{n_variants}\\n')
    handle.write(f'has_dosage\\t{"dosage" in mt.entry}\\n')
PY
    """
}


process QC_FILTER {
    tag "${cohort_id}"

    container params.hail_image
    cpus 8
    memory '32 GB'
    time '12h'

    input:
    path mt
    val cohort_id

    output:
    path 'qc.mt', emit: mt
    path 'sample-qc.tsv', emit: sample_qc
    path 'variant-qc.tsv', emit: variant_qc
    path 'qc-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    mkdir -p hail-tmp

    python3 <<'PY'
import hail as hl

hl.init(
    master='local[${task.cpus}]',
    tmp_dir='hail-tmp',
    local_tmpdir='hail-tmp',
    default_reference='${params.reference_genome}',
    spark_conf={'spark.driver.memory': '${params.hail_driver_memory}'},
    log='hail-qc.log',
)

mt = hl.read_matrix_table('${mt}')
variants_in, samples_in = mt.count()

# Entry-level filters only apply to called data carrying quality fields.
# Imputed BGEN dosages carry neither GQ nor DP, so these are skipped there.
entry_filters_applied = []
if 'GQ' in mt.entry and ${params.min_gq} > 0:
    mt = mt.filter_entries(hl.is_defined(mt.GQ) & (mt.GQ >= ${params.min_gq}))
    entry_filters_applied.append('GQ')
if 'DP' in mt.entry and ${params.min_dp} > 0:
    mt = mt.filter_entries(hl.is_defined(mt.DP) & (mt.DP >= ${params.min_dp}))
    entry_filters_applied.append('DP')

mt = hl.variant_qc(mt)
mt.rows().select('variant_qc').flatten().export('variant-qc.tsv')

mt = mt.filter_rows(
    (mt.variant_qc.call_rate >= ${params.min_call_rate_variant})
    & (hl.min(mt.variant_qc.AF) >= ${params.min_maf})
    & (mt.variant_qc.p_value_hwe >= ${params.min_hwe_p})
)
variants_after_variant_qc = mt.count_rows()
if variants_after_variant_qc == 0:
    raise SystemExit('variant QC removed every variant')

mt = hl.sample_qc(mt)
mt.cols().select('sample_qc').flatten().export('sample-qc.tsv')

mt = mt.filter_cols(mt.sample_qc.call_rate >= ${params.min_call_rate_sample})
samples_retained = mt.count_cols()
if samples_retained == 0:
    raise SystemExit('sample QC removed every sample')

# Allele frequencies shift once samples are dropped, so re-derive them.
mt = mt.drop('variant_qc')
mt = hl.variant_qc(mt)
mt = mt.filter_rows(hl.min(mt.variant_qc.AF) >= ${params.min_maf})
variants_retained = mt.count_rows()

if samples_retained < ${params.min_samples_retained}:
    raise SystemExit(
        f'only {samples_retained} samples survived QC, '
        'below min_samples_retained=${params.min_samples_retained}'
    )
if variants_retained < ${params.min_variants_retained}:
    raise SystemExit(
        f'only {variants_retained} variants survived QC, '
        'below min_variants_retained=${params.min_variants_retained}'
    )

mt.write('qc.mt', overwrite=True)

with open('qc-metrics.tsv', 'w') as handle:
    handle.write(f'entry_filters_applied\\t{",".join(entry_filters_applied) or "none"}\\n')
    handle.write(f'samples_in\\t{samples_in}\\n')
    handle.write(f'samples_retained\\t{samples_retained}\\n')
    handle.write(f'samples_removed\\t{samples_in - samples_retained}\\n')
    handle.write(f'variants_in\\t{variants_in}\\n')
    handle.write(f'variants_after_variant_qc\\t{variants_after_variant_qc}\\n')
    handle.write(f'variants_retained\\t{variants_retained}\\n')
    handle.write(f'variants_removed\\t{variants_in - variants_retained}\\n')
PY
    """
}


process COMPUTE_PCA {
    tag "${cohort_id}"

    container params.hail_image
    cpus 8
    memory '32 GB'
    time '8h'

    input:
    path mt
    val cohort_id

    output:
    path 'pcs.tsv', emit: pcs
    path 'pca-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    mkdir -p hail-tmp

    python3 <<'PY'
import hail as hl

hl.init(
    master='local[${task.cpus}]',
    tmp_dir='hail-tmp',
    local_tmpdir='hail-tmp',
    default_reference='${params.reference_genome}',
    spark_conf={'spark.driver.memory': '${params.hail_driver_memory}'},
    log='hail-pca.log',
)

mt = hl.read_matrix_table('${mt}')

# Common, biallelic, autosomal variants only: the usual basis for structure axes.
common = mt.filter_rows(
    (hl.min(mt.variant_qc.AF) > 0.05)
    & (hl.len(mt.alleles) == 2)
    & (mt.locus.in_autosome())
)
n_candidate = common.count_rows()
if n_candidate == 0:
    raise SystemExit('no common autosomal biallelic variants available for PCA')

if ${params.pca_max_variants} > 0 and n_candidate > ${params.pca_max_variants}:
    common = common.sample_rows(${params.pca_max_variants} / n_candidate, seed=0)

if ${params.ld_prune ? 'True' : 'False'}:
    pruned = hl.ld_prune(
        common.GT,
        r2=${params.ld_prune_r2},
        bp_window_size=${params.ld_prune_window},
    )
    common = common.filter_rows(hl.is_defined(pruned[common.row_key]))

n_pca_variants = common.count_rows()
if n_pca_variants < ${params.n_pcs} * 10:
    raise SystemExit(
        f'{n_pca_variants} variants remain for PCA, too few for ${params.n_pcs} components'
    )

eigenvalues, scores, _ = hl.hwe_normalized_pca(common.GT, k=${params.n_pcs})

scores = scores.annotate(
    **{f'PC{i + 1}': scores.scores[i] for i in range(${params.n_pcs})}
).drop('scores')
scores.export('pcs.tsv')

total = sum(eigenvalues)
with open('pca-metrics.tsv', 'w') as handle:
    handle.write(f'pca_candidate_variants\\t{n_candidate}\\n')
    handle.write(f'pca_variants_used\\t{n_pca_variants}\\n')
    handle.write(f'pca_components\\t${params.n_pcs}\\n')
    for i, value in enumerate(eigenvalues):
        handle.write(f'pc{i + 1}_variance_fraction\\t{value / total:.6f}\\n')
PY
    """
}


process RUN_ASSOCIATION {
    tag "${cohort_id}"

    container params.hail_image
    cpus 8
    memory '32 GB'
    time '12h'

    input:
    path mt
    path pcs
    path phenotypes
    val cohort_id
    val phenotype_sha256

    output:
    path 'sumstats.tsv', emit: sumstats
    path 'assoc-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    if [[ -n '${phenotype_sha256}' ]]; then
      echo '${phenotype_sha256}  ${phenotypes}' | sha256sum -c -
    fi

    mkdir -p hail-tmp

    python3 <<'PY'
import hail as hl
import numpy as np
from scipy.stats import chi2

hl.init(
    master='local[${task.cpus}]',
    tmp_dir='hail-tmp',
    local_tmpdir='hail-tmp',
    default_reference='${params.reference_genome}',
    spark_conf={'spark.driver.memory': '${params.hail_driver_memory}'},
    log='hail-assoc.log',
)

mt = hl.read_matrix_table('${mt}')

phenotypes = hl.import_table(
    '${phenotypes}',
    impute=True,
    key='${params.sample_id_column}',
    types={'${params.sample_id_column}': hl.tstr},
)
pcs = hl.import_table('${pcs}', impute=True, key='s', types={'s': hl.tstr})

mt = mt.annotate_cols(pheno=phenotypes[mt.s], pcs=pcs[mt.s])

samples_total = mt.count_cols()
mt = mt.filter_cols(hl.is_defined(mt.pheno) & hl.is_defined(mt.pcs))
samples_matched = mt.count_cols()
if samples_matched == 0:
    raise SystemExit(
        'no genotype sample ID matched the phenotype table; '
        'check sample_id_column and cohort membership'
    )

response = mt.pheno['${params.phenotype_column}']
samples_with_phenotype = mt.aggregate_cols(hl.agg.count_where(hl.is_defined(response)))
if samples_with_phenotype < ${params.min_samples_retained}:
    raise SystemExit(
        f'only {samples_with_phenotype} samples carry a non-missing phenotype, '
        'below min_samples_retained=${params.min_samples_retained}'
    )

covariate_names = [c for c in '${params.covariate_columns}'.split(',') if c]
covariates = [1.0]
covariates += [hl.float64(mt.pheno[name]) for name in covariate_names]
covariates += [mt.pcs[f'PC{i + 1}'] for i in range(${params.n_pcs})]

is_binary = '${params.phenotype_type}' == 'binary'
if is_binary:
    distinct = mt.aggregate_cols(hl.agg.collect_as_set(response))
    observed = {v for v in distinct if v is not None}
    if not observed <= {0, 1, True, False}:
        raise SystemExit(
            f'binary phenotype must be coded 0/1, observed values: {sorted(map(str, observed))}'
        )
    cases = mt.aggregate_cols(hl.agg.count_where(hl.bool(response)))
    controls = samples_with_phenotype - cases
    if cases == 0 or controls == 0:
        raise SystemExit(f'binary phenotype has {cases} cases and {controls} controls')
    result = hl.logistic_regression_rows(
        test='wald',
        y=hl.bool(response),
        x=mt.GT.n_alt_alleles(),
        covariates=covariates,
    )
else:
    cases = controls = 'NA'
    result = hl.linear_regression_rows(
        y=hl.float64(response),
        x=mt.GT.n_alt_alleles(),
        covariates=covariates,
    )

result = result.annotate(
    chrom=result.locus.contig,
    pos=result.locus.position,
    ref=result.alleles[0],
    alt=result.alleles[1],
)
result = result.key_by()
result = result.select(
    'chrom', 'pos', 'ref', 'alt', 'beta', 'standard_error', 'p_value'
)
result.export('sumstats.tsv')

pvalues = np.loadtxt('sumstats.tsv', delimiter='\\t', skiprows=1, usecols=6, ndmin=1)
pvalues = pvalues[np.isfinite(pvalues) & (pvalues > 0)]
if pvalues.size == 0:
    raise SystemExit('association produced no usable p-values')

lambda_gc = float(np.median(chi2.isf(pvalues, 1)) / chi2.ppf(0.5, 1))
genome_wide = int((pvalues < 5e-8).sum())
suggestive = int((pvalues < 1e-5).sum())

with open('assoc-metrics.tsv', 'w') as handle:
    handle.write(f'phenotype\\t${params.phenotype_column}\\n')
    handle.write(f'phenotype_type\\t${params.phenotype_type}\\n')
    handle.write(f'covariates\\t{",".join(covariate_names) or "none"}\\n')
    handle.write(f'pcs_used\\t${params.n_pcs}\\n')
    handle.write(f'samples_total\\t{samples_total}\\n')
    handle.write(f'samples_matched\\t{samples_matched}\\n')
    handle.write(f'samples_with_phenotype\\t{samples_with_phenotype}\\n')
    handle.write(f'cases\\t{cases}\\n')
    handle.write(f'controls\\t{controls}\\n')
    handle.write(f'variants_tested\\t{pvalues.size}\\n')
    handle.write(f'lambda_gc\\t{lambda_gc:.4f}\\n')
    handle.write(f'genome_wide_significant\\t{genome_wide}\\n')
    handle.write(f'suggestive\\t{suggestive}\\n')

if not (${params.min_lambda_gc} <= lambda_gc <= ${params.max_lambda_gc}):
    raise SystemExit(
        f'lambda_gc {lambda_gc:.4f} outside '
        '[${params.min_lambda_gc}, ${params.max_lambda_gc}]; '
        'residual structure or QC failure is likely'
    )
PY
    """
}


process PLOT_RESULTS {
    tag "${cohort_id}"

    container params.hail_image
    cpus 2
    memory '8 GB'
    time '1h'

    input:
    path sumstats
    val cohort_id

    output:
    path 'manhattan.html', emit: manhattan
    path 'qq.html', emit: qq
    path 'top-hits.tsv', emit: top_hits

    script:
    """
    set -euo pipefail

    python3 <<'PY'
import numpy as np
import pandas as pd
from bokeh.io import save, output_file
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure

frame = pd.read_csv('${sumstats}', sep='\\t')
frame = frame[np.isfinite(frame['p_value']) & (frame['p_value'] > 0)].copy()
frame['neglog10p'] = -np.log10(frame['p_value'])

frame.sort_values('p_value').head(100).to_csv('top-hits.tsv', sep='\\t', index=False)

# Manhattan: lay chromosomes end to end on a single cumulative axis.
def contig_sort_key(name):
    bare = str(name).replace('chr', '')
    order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
    return order.get(bare, int(bare) if bare.isdigit() else 99)

contigs = sorted(frame['chrom'].unique(), key=contig_sort_key)
offset = 0
offsets = {}
ticks = []
for contig in contigs:
    offsets[contig] = offset
    span = frame.loc[frame['chrom'] == contig, 'pos'].max()
    ticks.append((offset + span / 2, str(contig).replace('chr', '')))
    offset += span

frame['cumulative_pos'] = frame['pos'] + frame['chrom'].map(offsets)
frame['color'] = [
    '#3b6ea5' if contigs.index(c) % 2 == 0 else '#9fb8d4' for c in frame['chrom']
]

output_file('manhattan.html', title='${cohort_id} Manhattan')
manhattan = figure(
    width=1400,
    height=500,
    title='${cohort_id} — association',
    x_axis_label='chromosome',
    y_axis_label='-log10(p)',
    tools='pan,box_zoom,wheel_zoom,reset,save',
)
manhattan.scatter(
    x='cumulative_pos',
    y='neglog10p',
    color='color',
    size=4,
    source=ColumnDataSource(frame),
)
manhattan.line(
    x=[0, offset], y=[-np.log10(5e-8)] * 2, color='#c0392b', line_dash='dashed'
)
manhattan.xaxis.ticker = [position for position, _ in ticks]
manhattan.xaxis.major_label_overrides = {
    position: label for position, label in ticks
}
manhattan.xgrid.grid_line_color = None
save(manhattan)

# QQ: observed against the uniform expectation under the null.
observed = np.sort(frame['neglog10p'].values)[::-1]
expected = -np.log10(np.arange(1, observed.size + 1) / (observed.size + 1))

output_file('qq.html', title='${cohort_id} QQ')
qq = figure(
    width=600,
    height=600,
    title='${cohort_id} — QQ',
    x_axis_label='expected -log10(p)',
    y_axis_label='observed -log10(p)',
    tools='pan,box_zoom,wheel_zoom,reset,save',
)
qq.scatter(x=expected, y=observed, size=4, color='#3b6ea5')
limit = float(max(expected.max(), observed.max()))
qq.line(x=[0, limit], y=[0, limit], color='#c0392b', line_dash='dashed')
save(qq)
PY
    """
}


process BUILD_PROVENANCE {
    tag "${cohort_id}"

    container params.hail_image
    cpus 1
    memory '2 GB'
    time '30m'

    input:
    path import_metrics
    path qc_metrics
    path pca_metrics
    path assoc_metrics
    val cohort_id

    output:
    path "${cohort_id}.provenance.txt", emit: provenance

    script:
    """
    set -euo pipefail

    metric() {
      awk -F '\\t' -v key="\$1" '\$1 == key { print \$2; found=1; exit } END { if (!found) exit 1 }' "\$2"
    }

    cat > '${cohort_id}.provenance.txt' <<EOF
cohort              ${cohort_id}
date                \$(date -Iseconds)
reference_genome    ${params.reference_genome}
genotype_format     \$(metric genotype_format '${import_metrics}')
configured_image    ${params.hail_image}

samples_imported    \$(metric samples_imported '${import_metrics}')
variants_imported   \$(metric variants_imported '${import_metrics}')
entry_filters       \$(metric entry_filters_applied '${qc_metrics}')
samples_retained    \$(metric samples_retained '${qc_metrics}')
variants_retained   \$(metric variants_retained '${qc_metrics}')

qc_thresholds       min_gq=${params.min_gq} min_dp=${params.min_dp}
                    variant_call_rate=${params.min_call_rate_variant}
                    sample_call_rate=${params.min_call_rate_sample}
                    maf=${params.min_maf} hwe_p=${params.min_hwe_p}

pca_variants_used   \$(metric pca_variants_used '${pca_metrics}')
pcs_used            \$(metric pcs_used '${assoc_metrics}')

phenotype           \$(metric phenotype '${assoc_metrics}')
phenotype_type      \$(metric phenotype_type '${assoc_metrics}')
covariates          \$(metric covariates '${assoc_metrics}')
samples_analysed    \$(metric samples_with_phenotype '${assoc_metrics}')
cases               \$(metric cases '${assoc_metrics}')
controls            \$(metric controls '${assoc_metrics}')
variants_tested     \$(metric variants_tested '${assoc_metrics}')
lambda_gc           \$(metric lambda_gc '${assoc_metrics}')
genome_wide_hits    \$(metric genome_wide_significant '${assoc_metrics}')
suggestive_hits     \$(metric suggestive '${assoc_metrics}')

NOT VALIDATED FOR CLINICAL USE.
Missing: relatedness estimation (no kinship filter or mixed model), ancestry
projection onto an external reference panel, sex-check, per-chromosome X/Y
handling, and replication against external summary statistics.
A cohort of this size is underpowered for common-variant discovery; treat
results as replication or lookup, not discovery.
EOF
    """
}


process COLLECT_RESULTS {
    tag "${cohort_id}"

    container params.hail_image
    cpus 1
    memory '2 GB'
    time '30m'

    publishDir params.output_dir, mode: 'copy', overwrite: false

    input:
    path sumstats
    path top_hits
    path manhattan
    path qq
    path pcs
    path sample_qc
    path variant_qc
    path provenance
    path import_metrics
    path qc_metrics
    path pca_metrics
    path assoc_metrics
    val cohort_id

    output:
    path 'results/*', emit: files

    script:
    """
    set -euo pipefail

    mkdir -p results
    cp '${sumstats}' "results/${cohort_id}.sumstats.tsv"
    cp '${top_hits}' "results/${cohort_id}.top-hits.tsv"
    cp '${manhattan}' "results/${cohort_id}.manhattan.html"
    cp '${qq}' "results/${cohort_id}.qq.html"
    cp '${pcs}' "results/${cohort_id}.pcs.tsv"
    cp '${sample_qc}' "results/${cohort_id}.sample-qc.tsv"
    cp '${variant_qc}' "results/${cohort_id}.variant-qc.tsv"
    cp '${provenance}' "results/${cohort_id}.provenance.txt"
    cp '${import_metrics}' "results/${cohort_id}.import-metrics.tsv"
    cp '${qc_metrics}' "results/${cohort_id}.qc-metrics.tsv"
    cp '${pca_metrics}' "results/${cohort_id}.pca-metrics.tsv"
    cp '${assoc_metrics}' "results/${cohort_id}.assoc-metrics.tsv"
    """
}


workflow {
    if (!params.genotype_uri) {
        error 'genotype_uri is required'
    }
    if (!params.phenotype_uri) {
        error 'phenotype_uri is required'
    }
    if (!params.cohort_id) {
        error 'cohort_id is required'
    }
    if (!(params.cohort_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) {
        error 'cohort_id may contain only letters, numbers, period, underscore, and hyphen'
    }
    if (!params.output_dir) {
        error 'output_dir is required'
    }
    if (!(params.genotype_format in ['mt', 'vcf', 'bgen'])) {
        error 'genotype_format must be mt, vcf, or bgen'
    }
    if (!(params.reference_genome in ['GRCh37', 'GRCh38'])) {
        error 'reference_genome must be GRCh37 or GRCh38'
    }
    if (!params.phenotype_column) {
        error 'phenotype_column is required'
    }
    if (!(params.phenotype_type in ['quantitative', 'binary'])) {
        error 'phenotype_type must be quantitative or binary'
    }
    if (params.genotype_sha256 && !(params.genotype_sha256 ==~ /[0-9a-fA-F]{64}/)) {
        error 'genotype_sha256 must contain 64 hexadecimal characters'
    }
    if (params.phenotype_sha256 && !(params.phenotype_sha256 ==~ /[0-9a-fA-F]{64}/)) {
        error 'phenotype_sha256 must contain 64 hexadecimal characters'
    }

    // Column names reach Hail as Python identifiers; keep them inert.
    column_pattern = /[A-Za-z_][A-Za-z0-9_.]{0,63}/
    if (!(params.sample_id_column ==~ column_pattern)) {
        error 'sample_id_column must be a plain column name'
    }
    if (!(params.phenotype_column ==~ column_pattern)) {
        error 'phenotype_column must be a plain column name'
    }
    params.covariate_columns.toString().tokenize(',').each { name ->
        if (!(name ==~ column_pattern)) {
            error "covariate column '${name}' must be a plain column name"
        }
    }

    n_pcs = params.n_pcs as Integer
    if (n_pcs < 1 || n_pcs > 50) {
        error 'n_pcs must be between 1 and 50'
    }
    [
        'min_call_rate_variant': params.min_call_rate_variant,
        'min_call_rate_sample': params.min_call_rate_sample,
        'min_maf': params.min_maf,
        'min_hwe_p': params.min_hwe_p,
        'ld_prune_r2': params.ld_prune_r2,
    ].each { name, value ->
        def fraction = value as BigDecimal
        if (fraction < 0 || fraction > 1) {
            error "${name} must be between 0 and 1"
        }
    }
    if ((params.min_lambda_gc as BigDecimal) > (params.max_lambda_gc as BigDecimal)) {
        error 'min_lambda_gc must not exceed max_lambda_gc'
    }

    genotype = channel.fromPath(params.genotype_uri, checkIfExists: true)
    phenotypes = channel.fromPath(params.phenotype_uri, checkIfExists: true)
    bgen_sample = params.bgen_sample_uri
        ? channel.fromPath(params.bgen_sample_uri, checkIfExists: true)
        : channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)

    cohort_id = channel.value(params.cohort_id)
    genotype_sha256 = channel.value(params.genotype_sha256.toString().toLowerCase())
    phenotype_sha256 = channel.value(params.phenotype_sha256.toString().toLowerCase())

    IMPORT_GENOTYPES(genotype, bgen_sample, cohort_id, genotype_sha256)

    QC_FILTER(IMPORT_GENOTYPES.out.mt, cohort_id)

    COMPUTE_PCA(QC_FILTER.out.mt, cohort_id)

    RUN_ASSOCIATION(
        QC_FILTER.out.mt,
        COMPUTE_PCA.out.pcs,
        phenotypes,
        cohort_id,
        phenotype_sha256,
    )

    PLOT_RESULTS(RUN_ASSOCIATION.out.sumstats, cohort_id)

    BUILD_PROVENANCE(
        IMPORT_GENOTYPES.out.metrics,
        QC_FILTER.out.metrics,
        COMPUTE_PCA.out.metrics,
        RUN_ASSOCIATION.out.metrics,
        cohort_id,
    )

    COLLECT_RESULTS(
        RUN_ASSOCIATION.out.sumstats,
        PLOT_RESULTS.out.top_hits,
        PLOT_RESULTS.out.manhattan,
        PLOT_RESULTS.out.qq,
        COMPUTE_PCA.out.pcs,
        QC_FILTER.out.sample_qc,
        QC_FILTER.out.variant_qc,
        BUILD_PROVENANCE.out.provenance,
        IMPORT_GENOTYPES.out.metrics,
        QC_FILTER.out.metrics,
        COMPUTE_PCA.out.metrics,
        RUN_ASSOCIATION.out.metrics,
        cohort_id,
    )
}
