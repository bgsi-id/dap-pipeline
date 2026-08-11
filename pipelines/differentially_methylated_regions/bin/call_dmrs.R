#!/usr/bin/env Rscript

# Cohort DMR calling with DSS. The data path follows nf-core/methylong's
# population-scale DSS implementation: one chr/pos/N/X table per biological
# replicate, DMLtest, callDML, then callDMR.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13) {
  stop(paste(
    'expected matrix_dir, output_dir, partition, region, p_threshold,',
    'delta, min_length, min_cpgs, merge_distance, significant_fraction,',
    'equal_dispersion, smoothing, and ncores'
  ))
}

matrix_dir <- args[[1]]
out <- args[[2]]
partition <- args[[3]]
genomic_region <- args[[4]]
p_threshold <- as.numeric(args[[5]])
delta <- as.numeric(args[[6]])
min_length <- as.integer(args[[7]])
min_cpgs <- as.integer(args[[8]])
merge_distance <- as.integer(args[[9]])
significant_fraction <- as.numeric(args[[10]])
equal_dispersion <- tolower(args[[11]]) == 'true'
smoothing <- tolower(args[[12]]) == 'true'
ncores <- as.integer(args[[13]])

dir.create(out)

empty_outputs <- function(status, sites) {
  write.table(
    data.frame(
      chr = character(), start = integer(), end = integer(), length = integer(),
      nCG = integer(), meanMethy1 = numeric(), meanMethy2 = numeric(),
      diff.Methy = numeric(), areaStat = numeric()
    ),
    file.path(out, 'dmr-results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      chr = character(), pos = integer(), mu1 = numeric(), mu2 = numeric(),
      diff = numeric(), diff.se = numeric(), stat = numeric(), phi1 = numeric(),
      phi2 = numeric(), pval = numeric(), fdr = numeric()
    ),
    file.path(out, 'cpg-results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      chr = character(), pos = integer(), mu1 = numeric(), mu2 = numeric(),
      diff = numeric(), diff.se = numeric(), stat = numeric(), phi1 = numeric(),
      phi2 = numeric(), pval = numeric(), fdr = numeric()
    ),
    file.path(out, 'dml-results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(metric = c('partition', 'status', 'cpgs_available'), value = c(partition, status, sites)),
    file.path(out, 'dmr-metrics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(
      partition = character(), chrom = character(), position = integer(),
      group = character(), mean_beta = numeric()
    ),
    file.path(out, 'top-dmr-profile.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
  )
  png(file.path(out, 'top-dmr-profile.png'), width = 1400, height = 800, res = 140)
  plot.new()
  text(.5, .55, paste('No', partition, 'DMR plot available'), cex = 1.2)
  text(.5, .45, paste('Status:', status), cex = .9, col = '#666666')
  dev.off()
  file.copy(
    file.path(out, 'top-dmr-profile.png'), file.path(out, 'top-dmr-dss.png'),
    overwrite = TRUE
  )
}

required <- c('DSS', 'bsseq')
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) stop('DMR image is missing packages: ', paste(missing, collapse = ', '))
suppressPackageStartupMessages(library(parallel))

read_count_matrix <- function(name) {
  frame <- read.delim(gzfile(file.path(matrix_dir, name)), check.names = FALSE)
  if (!'probe_id' %in% names(frame)) stop(name, ' lacks probe_id')
  frame
}

coverage_frame <- read_count_matrix('coverage.tsv.gz')
modified_frame <- read_count_matrix('modified.tsv.gz')
probes <- read.delim(gzfile(file.path(matrix_dir, 'probes.tsv.gz')), check.names = FALSE)
metadata <- read.delim(
  file.path(matrix_dir, 'metadata.tsv'), check.names = FALSE, stringsAsFactors = FALSE
)

if (!identical(coverage_frame$probe_id, modified_frame$probe_id) ||
    !identical(coverage_frame$probe_id, probes$probe_id)) {
  stop('coverage, modified-count, and probe-coordinate tables disagree')
}
if (!identical(names(coverage_frame)[-1], names(modified_frame)[-1]) ||
    !identical(names(coverage_frame)[-1], metadata$sample_id)) {
  stop('count matrices and cohort membership disagree')
}
if (anyDuplicated(paste(probes$chrom, probes$position, sep = ':'))) {
  stop('DSS requires one CpG per genomic coordinate; duplicate coordinates remain after preprocessing')
}

groups <- factor(metadata$group, levels = c('control', 'patient'))
if (anyNA(groups) || any(table(groups) < 3)) {
  stop('control and patient each require at least three biological replicates')
}

coverage <- as.matrix(coverage_frame[, -1, drop = FALSE])
modified <- as.matrix(modified_frame[, -1, drop = FALSE])
storage.mode(coverage) <- 'double'
storage.mode(modified) <- 'double'
if (any(coverage <= 0, na.rm = TRUE) || any(modified < 0, na.rm = TRUE) ||
    any(modified > coverage, na.rm = TRUE)) {
  stop('invalid DSS counts: require N > 0 and 0 <= X <= N')
}

if (nrow(probes) < min_cpgs) {
  empty_outputs('insufficient_sites', nrow(probes))
  writeLines(c(
    paste('partition', partition, sep = '\t'),
    paste('genomic_region', ifelse(nzchar(genomic_region), genomic_region, 'genome-wide'), sep = '\t'),
    'method\tDSS DMLtest + callDML + callDMR',
    'status\tinsufficient_sites',
    'not_clinical_use\ttrue'
  ), file.path(out, 'provenance.tsv'))
  quit(save = 'no', status = 0)
}

sample_names <- metadata$sample_id
sample_tables <- lapply(seq_along(sample_names), function(index) {
  keep <- is.finite(coverage[, index]) & is.finite(modified[, index])
  data.frame(
    chr = as.character(probes$chrom[keep]),
    pos = as.integer(probes$position[keep]),
    N = as.integer(coverage[keep, index]),
    X = as.integer(modified[keep, index])
  )
})
if (any(vapply(sample_tables, nrow, integer(1)) == 0)) {
  stop('at least one biological replicate has no retained CpGs')
}

bs_object <- DSS::makeBSseqData(sample_tables, sample_names)
patient_samples <- sample_names[groups == 'patient']
control_samples <- sample_names[groups == 'control']

# Group 1 is patient and group 2 is control, matching methylong's case/control
# ordering. Therefore DSS diff and DMR diff.Methy are patient minus control.
dml_test <- DSS::DMLtest(
  bs_object,
  group1 = patient_samples,
  group2 = control_samples,
  equal.disp = equal_dispersion,
  smoothing = smoothing,
  ncores = ncores
)
dml_calls <- DSS::callDML(dml_test, delta = delta, p.threshold = p_threshold)
dmr_calls <- DSS::callDMR(
  dml_test,
  delta = delta,
  p.threshold = p_threshold,
  minlen = min_length,
  minCG = min_cpgs,
  dis.merge = merge_distance,
  pct.sig = significant_fraction
)

write.table(
  dml_test, file.path(out, 'cpg-results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
)
write.table(
  dml_calls, file.path(out, 'dml-results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
)
write.table(
  dmr_calls, file.path(out, 'dmr-results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
)

metrics <- data.frame(
  metric = c(
    'partition', 'method', 'samples', 'control_samples', 'patient_samples',
    'cpgs_tested', 'dmls', 'dmrs'
  ),
  value = c(
    partition, 'DSS', length(sample_names), length(control_samples), length(patient_samples),
    nrow(dml_test), nrow(dml_calls), nrow(dmr_calls)
  )
)
write.table(
  metrics, file.path(out, 'dmr-metrics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
)

profile <- data.frame(
  partition = character(), chrom = character(), position = integer(),
  group = character(), mean_beta = numeric()
)
if (nrow(dmr_calls)) {
  png(file.path(out, 'top-dmr-dss.png'), width = 1400, height = 800, res = 140)
  DSS::showOneDMR(dmr_calls[1, , drop = FALSE], bs_object)
  dev.off()

  top <- dmr_calls[1, , drop = FALSE]
  in_top <- probes$chrom == top$chr[[1]] &
    probes$position >= top$start[[1]] & probes$position <= top$end[[1]]
  beta <- modified / coverage
  group_profiles <- list()
  for (group_name in c('control', 'patient')) {
    members <- groups == group_name
    group_profile <- data.frame(
      partition = partition,
      chrom = as.character(probes$chrom[in_top]),
      position = as.integer(probes$position[in_top]),
      group = group_name,
      mean_beta = rowMeans(beta[in_top, members, drop = FALSE], na.rm = TRUE)
    )
    profile <- rbind(profile, group_profile)
    group_profiles[[group_name]] <- group_profile
  }

  png(file.path(out, 'top-dmr-profile.png'), width = 1400, height = 800, res = 140)
  colours <- c(control = '#0072B2', patient = '#D55E00')
  positions <- sort(unique(profile$position))
  plot(
    range(positions), c(0, 1), type = 'n', las = 1,
    xlab = paste0(top$chr[[1]], ' genomic position (GRCh38)'),
    ylab = 'Mean methylation fraction',
    main = paste0('Top DSS DMR — ', partition, ' cohort profiles')
  )
  grid(col = '#E5E5E5')
  for (group_name in names(group_profiles)) {
    item <- group_profiles[[group_name]]
    points(item$position, item$mean_beta, pch = 16, cex = .65,
           col = adjustcolor(colours[[group_name]], alpha.f = .4))
    if (nrow(item) >= 4 && length(unique(item$position)) >= 4) {
      fitted <- smooth.spline(item$position, item$mean_beta, spar = 0.6)
      curve <- predict(fitted, seq(min(item$position), max(item$position), length.out = 300))
      lines(curve$x, pmin(1, pmax(0, curve$y)), col = colours[[group_name]], lwd = 3)
    } else {
      lines(item$position, item$mean_beta, col = colours[[group_name]], lwd = 3)
    }
  }
  legend(
    'topright', legend = c(
      paste0('control (n=', length(control_samples), ')'),
      paste0('patient (n=', length(patient_samples), ')')
    ), col = colours, lwd = 3, bty = 'n'
  )
  mtext(
    sprintf('DSS patient - control = %.3f; curves are display-only cohort means.', top$diff.Methy[[1]]),
    side = 1, line = 4, cex = .8
  )
  dev.off()
} else {
  png(file.path(out, 'top-dmr-profile.png'), width = 1400, height = 800, res = 140)
  plot.new()
  text(.5, .55, paste('No', partition, 'DMR passed DSS thresholds'), cex = 1.2)
  text(.5, .45, 'No regional profile is available', cex = .9, col = '#666666')
  dev.off()
  file.copy(
    file.path(out, 'top-dmr-profile.png'), file.path(out, 'top-dmr-dss.png'),
    overwrite = TRUE
  )
}
write.table(
  profile, file.path(out, 'top-dmr-profile.tsv'), sep = '\t', quote = FALSE, row.names = FALSE
)

writeLines(c(
  paste('partition', partition, sep = '\t'),
  paste('genomic_region', ifelse(nzchar(genomic_region), genomic_region, 'genome-wide'), sep = '\t'),
  'contrast\tpatient - control',
  'method\tDSS DMLtest + callDML + callDMR',
  paste('p_threshold', p_threshold, sep = '\t'),
  paste('delta', delta, sep = '\t'),
  paste('min_length', min_length, sep = '\t'),
  paste('min_cpgs', min_cpgs, sep = '\t'),
  paste('merge_distance', merge_distance, sep = '\t'),
  paste('significant_fraction', significant_fraction, sep = '\t'),
  paste('equal_dispersion', equal_dispersion, sep = '\t'),
  paste('smoothing', smoothing, sep = '\t'),
  'input\tmodkit bedMethyl valid coverage (N) and modified counts (X)',
  'reference_genome\tGRCh38',
  if (partition == 'combined') 'interpretation\tprimary all-read result' else 'interpretation\texploratory haplotype/unphased partition; HP labels are not aligned across samples',
  'not_clinical_use\ttrue'
), file.path(out, 'provenance.tsv'))
