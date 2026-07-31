nextflow.enable.dsl = 2

/*
 * Cohort polygenic score (Hail)
 *
 * import -> parse PGS Catalog weights -> allele-aware match -> weighted sum
 *
 * Scoring is a dosage-weighted sum over variants that matched the scoring file
 * on locus and alleles. Coverage is measured and gated: a score computed from a
 * fraction of its weights is not a smaller score, it is a different score, and
 * the run fails rather than publishing one.
 */

params.genotype_uri = null
params.genotype_format = 'mt'
params.genotype_sha256 = ''
params.bgen_sample_uri = ''
params.score_uri = null
params.score_sha256 = ''

params.cohort_id = null
params.score_id = null
params.output_dir = 'results'
params.reference_genome = 'GRCh38'

params.exclude_ambiguous = true
params.min_variant_match_rate = 0.9
params.min_samples_scored = 100
params.max_sample_missing_rate = 0.1

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

# Imputed dosages carry more information than hard calls and are preferred
# when both are present.
if 'dosage' in mt.entry:
    dosage_source = 'dosage'
elif 'GT' in mt.entry:
    dosage_source = 'GT'
else:
    raise SystemExit('input has neither dosage nor GT entry field')

mt = mt.key_rows_by('locus', 'alleles')

n_variants, n_samples = mt.count()
if n_samples == 0:
    raise SystemExit('input contains no samples')
if n_variants == 0:
    raise SystemExit('input contains no variants')

mt.write('cohort.mt', overwrite=True)

with open('import-metrics.tsv', 'w') as handle:
    handle.write(f'genotype_format\\t{fmt}\\n')
    handle.write(f'dosage_source\\t{dosage_source}\\n')
    handle.write(f'samples_imported\\t{n_samples}\\n')
    handle.write(f'variants_imported\\t{n_variants}\\n')
PY
    """
}


process IMPORT_SCORE {
    tag "${score_id}"

    container params.hail_image
    cpus 1
    memory '4 GB'
    time '1h'

    input:
    path score_file
    val score_id
    val score_sha256

    output:
    path 'weights.tsv', emit: weights
    path 'score-metrics.tsv', emit: metrics

    script:
    """
    set -euo pipefail

    if [[ -n '${score_sha256}' ]]; then
      echo '${score_sha256}  ${score_file}' | sha256sum -c -
    fi

    python3 <<'PY'
import gzip
import pandas as pd

path = '${score_file}'
opener = gzip.open if path.endswith('.gz') else open
with opener(path, 'rt') as handle:
    frame = pd.read_csv(handle, sep='\\t', comment='#', dtype=str)

# PGS Catalog harmonised files carry hm_* positions already lifted to the
# target build. Prefer them; fall back to the author-reported coordinates.
def pick(primary, fallback):
    if primary in frame.columns:
        series = frame[primary]
        if fallback in frame.columns:
            return series.fillna(frame[fallback])
        return series
    if fallback in frame.columns:
        return frame[fallback]
    return None

chrom = pick('hm_chr', 'chr_name')
position = pick('hm_pos', 'chr_position')
other = pick('other_allele', 'hm_inferOtherAllele')

missing = [
    name
    for name, series in (
        ('chromosome', chrom),
        ('position', position),
        ('other_allele', other),
    )
    if series is None
]
if missing:
    raise SystemExit(f'scoring file lacks required columns: {missing}')
if 'effect_allele' not in frame.columns or 'effect_weight' not in frame.columns:
    raise SystemExit('scoring file lacks effect_allele or effect_weight')

weights = pd.DataFrame(
    {
        'chrom': chrom,
        'pos': position,
        'effect_allele': frame['effect_allele'],
        'other_allele': other,
        'weight': frame['effect_weight'],
    }
)
rows_in = len(weights)

weights = weights.dropna()
weights = weights[weights['pos'].str.strip().str.isdigit()]
weights['weight'] = pd.to_numeric(weights['weight'], errors='coerce')
weights = weights.dropna(subset=['weight'])

for column in ('effect_allele', 'other_allele'):
    weights[column] = weights[column].str.strip().str.upper()
    weights = weights[weights[column].str.fullmatch('[ACGT]+')]

# Hail's built-in GRCh38 uses chr-prefixed contigs; GRCh37 does not.
weights['chrom'] = weights['chrom'].str.strip().str.replace('^chr', '', regex=True)
if '${params.reference_genome}' == 'GRCh38':
    weights['chrom'] = 'chr' + weights['chrom']

weights = weights.drop_duplicates(subset=['chrom', 'pos', 'effect_allele'])
rows_usable = len(weights)
if rows_usable == 0:
    raise SystemExit('scoring file produced no usable weights after parsing')

weights.to_csv('weights.tsv', sep='\\t', index=False)

with open('score-metrics.tsv', 'w') as handle:
    handle.write(f'score_id\\t${score_id}\\n')
    handle.write(f'weights_in_file\\t{rows_in}\\n')
    handle.write(f'weights_usable\\t{rows_usable}\\n')
    handle.write(f'weights_discarded\\t{rows_in - rows_usable}\\n')
PY
    """
}


process SCORE_SAMPLES {
    tag "${cohort_id}"

    container params.hail_image
    cpus 8
    memory '32 GB'
    time '8h'

    input:
    path mt
    path weights
    path score_metrics
    val cohort_id

    output:
    path 'scores.tsv', emit: scores
    path 'match-metrics.tsv', emit: metrics
    path 'unmatched-weights.tsv', emit: unmatched

    script:
    """
    set -euo pipefail

    mkdir -p hail-tmp

    python3 <<'PY'
import hail as hl
import pandas as pd

hl.init(
    master='local[${task.cpus}]',
    tmp_dir='hail-tmp',
    local_tmpdir='hail-tmp',
    default_reference='${params.reference_genome}',
    spark_conf={'spark.driver.memory': '${params.hail_driver_memory}'},
    log='hail-score.log',
)

mt = hl.read_matrix_table('${mt}')
mt = mt.filter_rows(hl.len(mt.alleles) == 2)

base = hl.import_table(
    '${weights}',
    types={
        'chrom': hl.tstr,
        'pos': hl.tint32,
        'effect_allele': hl.tstr,
        'other_allele': hl.tstr,
        'weight': hl.tfloat64,
    },
)
base = base.annotate(
    locus=hl.locus(base.chrom, base.pos, reference_genome='${params.reference_genome}')
)
weights_total = base.count()

# Palindromic sites cannot be oriented by allele identity alone. Excluding them
# is the safe default; including them silently risks sign-flipped weights.
palindromic = hl.literal({'AT', 'TA', 'CG', 'GC'})
base = base.annotate(
    ambiguous=palindromic.contains(base.effect_allele + base.other_allele)
)
weights_ambiguous = base.aggregate(hl.agg.count_where(base.ambiguous))
if ${params.exclude_ambiguous ? 'True' : 'False'}:
    base = base.filter(~base.ambiguous)

# Emit both allele orientations so a single key join resolves which of the
# cohort's ref/alt is the effect allele.
forward = base.annotate(alleles=[base.other_allele, base.effect_allele], flip=False)
reverse = base.annotate(alleles=[base.effect_allele, base.other_allele], flip=True)
oriented = forward.union(reverse).key_by('locus', 'alleles')

mt = mt.annotate_rows(w=oriented[mt.row_key])
mt = mt.filter_rows(hl.is_defined(mt.w))
variants_matched = mt.count_rows()

match_rate = variants_matched / weights_total if weights_total else 0.0

rows = mt.rows()
matched_frame = rows.select(effect_allele=rows.w.effect_allele).key_by().to_pandas()
matched_pairs = set(
    zip(
        matched_frame['locus'].astype(str),
        matched_frame['effect_allele'].astype(str),
    )
)
weight_frame = pd.read_csv('${weights}', sep='\\t')
weight_frame['locus_key'] = (
    weight_frame['chrom'].astype(str) + ':' + weight_frame['pos'].astype(str)
)
weight_frame['matched'] = [
    (locus, allele) in matched_pairs
    for locus, allele in zip(weight_frame['locus_key'], weight_frame['effect_allele'])
]
weight_frame[~weight_frame['matched']].drop(columns=['matched']).to_csv(
    'unmatched-weights.tsv', sep='\\t', index=False
)

with open('match-metrics.tsv', 'w') as handle:
    handle.write(f'weights_total\\t{weights_total}\\n')
    handle.write(f'weights_ambiguous\\t{weights_ambiguous}\\n')
    handle.write(f'exclude_ambiguous\\t${params.exclude_ambiguous}\\n')
    handle.write(f'variants_matched\\t{variants_matched}\\n')
    handle.write(f'variant_match_rate\\t{match_rate:.6f}\\n')

if match_rate < ${params.min_variant_match_rate}:
    raise SystemExit(
        f'variant match rate {match_rate:.4f} below '
        'min_variant_match_rate=${params.min_variant_match_rate}; '
        'the cohort does not carry enough of this score to compute it'
    )

if 'dosage' in mt.entry:
    alt_dosage = mt.dosage
else:
    alt_dosage = mt.GT.n_alt_alleles()

mt = mt.annotate_entries(
    effect_dosage=hl.if_else(mt.w.flip, 2 - alt_dosage, alt_dosage)
)

# Cohort allele frequency is the mean-imputation target for missing calls.
mt = mt.annotate_rows(effect_af=hl.agg.mean(mt.effect_dosage) / 2)

mt = mt.annotate_cols(
    raw_score=hl.agg.sum(
        mt.w.weight * hl.or_else(mt.effect_dosage, 2 * mt.effect_af)
    ),
    genotypes_missing=hl.agg.count_where(hl.is_missing(mt.effect_dosage)),
)

scored = mt.cols().select('raw_score', 'genotypes_missing')
scored = scored.annotate(variants_used=variants_matched)
scored.export('scores.tsv')
PY

    python3 <<'PY'
import pandas as pd

frame = pd.read_csv('scores.tsv', sep='\\t')
frame = frame.rename(columns={'s': 'sample_id'})

samples_scored = len(frame)
if samples_scored < ${params.min_samples_scored}:
    raise SystemExit(
        f'only {samples_scored} samples scored, '
        'below min_samples_scored=${params.min_samples_scored}'
    )

frame['sample_missing_rate'] = frame['genotypes_missing'] / frame['variants_used']
worst = frame['sample_missing_rate'].max()
if worst > ${params.max_sample_missing_rate}:
    offenders = int((frame['sample_missing_rate'] > ${params.max_sample_missing_rate}).sum())
    raise SystemExit(
        f'{offenders} sample(s) exceed max_sample_missing_rate='
        '${params.max_sample_missing_rate} (worst {worst:.4f}); '
        'their scores would be dominated by imputed dosages'
    )

# Cohort-internal standardisation only. See the provenance caveat: this is a
# rank within this cohort, not a population reference percentile.
mean = frame['raw_score'].mean()
sd = frame['raw_score'].std(ddof=1)
frame['z_score'] = (frame['raw_score'] - mean) / sd if sd and sd > 0 else float('nan')
frame['cohort_percentile'] = frame['raw_score'].rank(pct=True) * 100

frame.to_csv('scores.tsv', sep='\\t', index=False)

with open('match-metrics.tsv', 'a') as handle:
    handle.write(f'samples_scored\\t{samples_scored}\\n')
    handle.write(f'score_mean\\t{mean:.6f}\\n')
    handle.write(f'score_sd\\t{sd:.6f}\\n')
    handle.write(f'max_sample_missing_rate_observed\\t{worst:.6f}\\n')
PY
    """
}


process PLOT_SCORES {
    tag "${cohort_id}"

    container params.hail_image
    cpus 2
    memory '8 GB'
    time '1h'

    input:
    path scores
    val cohort_id

    output:
    path 'score-distribution.html', emit: distribution

    script:
    """
    set -euo pipefail

    python3 <<'PY'
import numpy as np
import pandas as pd
from bokeh.io import output_file, save
from bokeh.layouts import column
from bokeh.plotting import figure

frame = pd.read_csv('${scores}', sep='\\t')

output_file('score-distribution.html', title='${cohort_id} polygenic score')

counts, edges = np.histogram(frame['raw_score'].values, bins=50)
histogram = figure(
    width=900,
    height=400,
    title='${cohort_id} — raw score distribution',
    x_axis_label='raw polygenic score',
    y_axis_label='samples',
    tools='pan,box_zoom,wheel_zoom,reset,save',
)
histogram.quad(
    top=counts, bottom=0, left=edges[:-1], right=edges[1:],
    fill_color='#3b6ea5', line_color='white',
)

ordered = frame.sort_values('raw_score')
quantiles = figure(
    width=900,
    height=400,
    title='${cohort_id} — score by cohort percentile',
    x_axis_label='cohort percentile',
    y_axis_label='raw polygenic score',
    tools='pan,box_zoom,wheel_zoom,reset,save',
)
quantiles.scatter(
    x=ordered['cohort_percentile'], y=ordered['raw_score'], size=4, color='#3b6ea5'
)

save(column(histogram, quantiles))
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
    path score_metrics
    path match_metrics
    val cohort_id
    val score_id

    output:
    path "${cohort_id}.provenance.txt", emit: provenance

    script:
    """
    set -euo pipefail

    metric() {
      awk -F '\\t' -v key="\$1" '\$1 == key { print \$2; found=1; exit } END { if (!found) exit 1 }' "\$2"
    }

    cat > '${cohort_id}.provenance.txt' <<EOF
cohort               ${cohort_id}
score                ${score_id}
date                 \$(date -Iseconds)
reference_genome     ${params.reference_genome}
genotype_format      \$(metric genotype_format '${import_metrics}')
dosage_source        \$(metric dosage_source '${import_metrics}')
configured_image     ${params.hail_image}

samples_imported     \$(metric samples_imported '${import_metrics}')
variants_imported    \$(metric variants_imported '${import_metrics}')

weights_in_file      \$(metric weights_in_file '${score_metrics}')
weights_usable       \$(metric weights_usable '${score_metrics}')
weights_discarded    \$(metric weights_discarded '${score_metrics}')
weights_ambiguous    \$(metric weights_ambiguous '${match_metrics}')
exclude_ambiguous    \$(metric exclude_ambiguous '${match_metrics}')
variants_matched     \$(metric variants_matched '${match_metrics}')
variant_match_rate   \$(metric variant_match_rate '${match_metrics}')
min_match_rate       ${params.min_variant_match_rate}

samples_scored       \$(metric samples_scored '${match_metrics}')
score_mean           \$(metric score_mean '${match_metrics}')
score_sd             \$(metric score_sd '${match_metrics}')
worst_sample_missing \$(metric max_sample_missing_rate_observed '${match_metrics}')

NOT VALIDATED FOR CLINICAL USE.
cohort_percentile ranks each sample against this cohort only. It is not a
population reference percentile and is not ancestry-matched. Two cohorts scored
with the same weights are not comparable on this column, and a percentile from
a cohort whose ancestry composition differs from the score's discovery
population does not carry the score's published effect size.
Missing: ancestry-matched reference distribution, projection of cohort samples
onto reference PCs, per-ancestry calibration, and validation against a scored
control cohort.
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
    path scores
    path distribution
    path unmatched
    path provenance
    path import_metrics
    path score_metrics
    path match_metrics
    val cohort_id

    output:
    path 'results/*', emit: files

    script:
    """
    set -euo pipefail

    mkdir -p results
    cp '${scores}' "results/${cohort_id}.scores.tsv"
    cp '${distribution}' "results/${cohort_id}.score-distribution.html"
    cp '${unmatched}' "results/${cohort_id}.unmatched-weights.tsv"
    cp '${provenance}' "results/${cohort_id}.provenance.txt"
    cp '${import_metrics}' "results/${cohort_id}.import-metrics.tsv"
    cp '${score_metrics}' "results/${cohort_id}.score-metrics.tsv"
    cp '${match_metrics}' "results/${cohort_id}.match-metrics.tsv"
    """
}


workflow {
    if (!params.genotype_uri) {
        error 'genotype_uri is required'
    }
    if (!params.score_uri) {
        error 'score_uri is required'
    }
    if (!params.cohort_id) {
        error 'cohort_id is required'
    }
    if (!params.score_id) {
        error 'score_id is required to identify the weights in provenance'
    }
    safe_id = /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/
    if (!(params.cohort_id ==~ safe_id)) {
        error 'cohort_id may contain only letters, numbers, period, underscore, and hyphen'
    }
    if (!(params.score_id ==~ safe_id)) {
        error 'score_id may contain only letters, numbers, period, underscore, and hyphen'
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
    if (params.genotype_sha256 && !(params.genotype_sha256 ==~ /[0-9a-fA-F]{64}/)) {
        error 'genotype_sha256 must contain 64 hexadecimal characters'
    }
    if (params.score_sha256 && !(params.score_sha256 ==~ /[0-9a-fA-F]{64}/)) {
        error 'score_sha256 must contain 64 hexadecimal characters'
    }

    match_rate = params.min_variant_match_rate as BigDecimal
    if (match_rate < 0 || match_rate > 1) {
        error 'min_variant_match_rate must be between 0 and 1'
    }
    missing_rate = params.max_sample_missing_rate as BigDecimal
    if (missing_rate < 0 || missing_rate > 1) {
        error 'max_sample_missing_rate must be between 0 and 1'
    }

    genotype = channel.fromPath(params.genotype_uri, checkIfExists: true)
    score_file = channel.fromPath(params.score_uri, checkIfExists: true)
    bgen_sample = params.bgen_sample_uri
        ? channel.fromPath(params.bgen_sample_uri, checkIfExists: true)
        : channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)

    cohort_id = channel.value(params.cohort_id)
    score_id = channel.value(params.score_id)
    genotype_sha256 = channel.value(params.genotype_sha256.toString().toLowerCase())
    score_sha256 = channel.value(params.score_sha256.toString().toLowerCase())

    IMPORT_GENOTYPES(genotype, bgen_sample, cohort_id, genotype_sha256)

    IMPORT_SCORE(score_file, score_id, score_sha256)

    SCORE_SAMPLES(
        IMPORT_GENOTYPES.out.mt,
        IMPORT_SCORE.out.weights,
        IMPORT_SCORE.out.metrics,
        cohort_id,
    )

    PLOT_SCORES(SCORE_SAMPLES.out.scores, cohort_id)

    BUILD_PROVENANCE(
        IMPORT_GENOTYPES.out.metrics,
        IMPORT_SCORE.out.metrics,
        SCORE_SAMPLES.out.metrics,
        cohort_id,
        score_id,
    )

    COLLECT_RESULTS(
        SCORE_SAMPLES.out.scores,
        PLOT_SCORES.out.distribution,
        SCORE_SAMPLES.out.unmatched,
        BUILD_PROVENANCE.out.provenance,
        IMPORT_GENOTYPES.out.metrics,
        IMPORT_SCORE.out.metrics,
        SCORE_SAMPLES.out.metrics,
        cohort_id,
    )
}
