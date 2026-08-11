nextflow.enable.dsl = 2

include { WRITE_RESULTS } from '../base_nf/modules/bundle'

params.methylation_matrix_uri = null
params.phenotype_uri = null
params.probe_manifest_uri = null
params.cohort_id = null
params.dap_output_uri = null
params.sample_id_column = 'sample_id'
params.phenotype_column = null
params.covariate_columns = ''
params.fdr_threshold = 0.05
params.lambda = 1000
params.bandwidth_scaling = 2
params.min_cpgs = 2
params.min_delta_beta = 0.0
params.profile_spar = 0.6
params.dmr_memory = '32 GB'

process VALIDATE_DMR_INPUTS {
    tag params.cohort_id
    cpus 1
    memory '4 GB'
    container params.r_image

    input:
    path methylation
    path phenotypes
    path manifest

    output:
    path 'validated'

    script:
    """
    set -euo pipefail
    mkdir validated
    Rscript - '${methylation}' '${phenotypes}' '${manifest}' validated <<'RS'
args <- commandArgs(trailingOnly=TRUE); out <- args[[4]]
read_auto <- function(path) read.delim(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep=ifelse(grepl('[.]csv\$',path,ignore.case=TRUE),',','\t'))
beta <- read_auto(args[[1]]); pheno <- read_auto(args[[2]]); probes <- read_auto(args[[3]])
if (ncol(beta) < 3) stop('methylation matrix requires probe_id plus at least two samples')
names(beta)[1] <- 'probe_id'
required_probe <- c('probe_id','chrom','position')
if (!all(required_probe %in% names(probes))) stop('probe manifest requires probe_id, chrom, and position')
sid <- '${params.sample_id_column}'; trait <- '${params.phenotype_column}'
if (!all(c(sid,trait) %in% names(pheno))) stop('phenotype file is missing sample or phenotype column')
if (anyDuplicated(beta\$probe_id) || anyDuplicated(probes\$probe_id) || anyDuplicated(pheno[[sid]])) stop('probe and sample identifiers must be unique')
sample_ids <- intersect(names(beta)[-1], as.character(pheno[[sid]]))
if (length(sample_ids) < 20) stop('DMR analysis requires at least 20 matched samples')
pheno <- pheno[match(sample_ids,pheno[[sid]]),,drop=FALSE]
groups <- factor(pheno[[trait]])
if (nlevels(groups) != 2) stop('DMR analysis requires exactly two phenotype groups')
if (any(table(groups) < 3)) stop('DMR analysis requires at least three samples in each phenotype group')
values <- as.matrix(beta[,sample_ids,drop=FALSE]); storage.mode(values) <- 'double'
if (any(!is.finite(values),na.rm=TRUE) || any(values < 0 | values > 1,na.rm=TRUE)) stop('methylation values must be beta values in [0,1]')
keep <- rowMeans(is.na(values)) <= 0.05 & beta\$probe_id %in% probes\$probe_id
values <- values[keep,,drop=FALSE]; rownames(values) <- beta\$probe_id[keep]
if (nrow(values) < 100) stop('fewer than 100 located probes remain after filtering')
for (i in seq_len(nrow(values))) values[i,is.na(values[i,])] <- median(values[i,],na.rm=TRUE)
saveRDS(list(beta=values,pheno=pheno,groups=groups,sample_ids=sample_ids),file.path(out,'inputs.rds'))
write.table(probes[match(rownames(values),probes\$probe_id),,drop=FALSE],file.path(out,'probes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(metric=c('samples','probes','group_1','group_2'),value=c(length(sample_ids),nrow(values),sum(groups==levels(groups)[1]),sum(groups==levels(groups)[2]))),file.path(out,'validation-metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
RS
    """

    stub:
    """
    mkdir validated; touch validated/inputs.rds
    printf 'probe_id\tchrom\tposition\n' > validated/probes.tsv
    printf 'metric\tvalue\nsamples\t20\nprobes\t100\n' > validated/validation-metrics.tsv
    """
}

process CALL_DMRS {
    tag params.cohort_id
    cpus 4
    memory params.dmr_memory
    container params.r_image

    input:
    path validated

    output:
    path 'dmr'

    script:
    """
    set -euo pipefail
    mkdir dmr
    Rscript - validated dmr <<'RS'
args <- commandArgs(trailingOnly=TRUE); source <- args[[1]]; out <- args[[2]]
required <- c('limma','DMRcate','GenomicRanges','IRanges')
missing <- required[!vapply(required,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))]
if (length(missing)) stop('DMR image is missing packages: ',paste(missing,collapse=', '))
data <- readRDS(file.path(source,'inputs.rds')); probes <- read.delim(file.path(source,'probes.tsv'),check.names=FALSE)
beta <- data\$beta; groups <- data\$groups; pheno <- data\$pheno
covariates <- Filter(nzchar,strsplit('${params.covariate_columns}',',',fixed=TRUE)[[1]])
if (!all(covariates %in% names(pheno))) stop('covariate column is missing')
design_data <- data.frame(group=groups,pheno[,covariates,drop=FALSE],check.names=FALSE)
design <- model.matrix(~.,data=design_data); coef_index <- grep('^group',colnames(design))[1]
if (is.na(coef_index)) stop('phenotype coefficient is absent from design')
clipped <- pmin(pmax(beta,1e-6),1-1e-6); mvalues <- log2(clipped/(1-clipped))
fit <- limma::eBayes(limma::lmFit(mvalues,design),robust=TRUE)
p <- fit\$p.value[,coef_index]; q <- p.adjust(p,'BH'); statistic <- fit\$t[,coef_index]
delta <- rowMeans(beta[,groups==levels(groups)[2],drop=FALSE])-rowMeans(beta[,groups==levels(groups)[1],drop=FALSE])
ranges <- GenomicRanges::GRanges(seqnames=probes\$chrom,ranges=IRanges::IRanges(as.integer(probes\$position),as.integer(probes\$position)),stat=statistic,rawpval=p,diff=delta,ind.fdr=q,is.sig=q<${params.fdr_threshold})
names(ranges) <- rownames(beta)
annotated <- methods::new('CpGannotated',ranges=ranges)
called <- DMRcate::dmrcate(annotated,lambda=${params.lambda},C=${params.bandwidth_scaling},min.cpgs=${params.min_cpgs},betacutoff=${params.min_delta_beta})
dmrs <- as.data.frame(DMRcate::extractRanges(called,genome='hg38'))
write.table(data.frame(probe_id=rownames(beta),chrom=probes\$chrom,position=probes\$position,delta_beta=delta,t=statistic,p_value=p,fdr=q),file.path(out,'cpg-results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(dmrs,file.path(out,'dmr-results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
metrics <- data.frame(metric=c('samples','cpgs_tested','cpgs_fdr','dmrs','group_reference','group_comparison'),value=c(ncol(beta),nrow(beta),sum(q<${params.fdr_threshold},na.rm=TRUE),nrow(dmrs),levels(groups)[1],levels(groups)[2]))
write.table(metrics,file.path(out,'dmr-metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

# Methylartist-inspired summary for the highest-ranked DMR: show only the
# normalized group profiles and uncertainty bands, never per-sample pileups.
profile_columns <- c('chrom','position','group','mean_beta','ci_lower','ci_upper','smoothed_beta')
profile <- data.frame(matrix(ncol=length(profile_columns),nrow=0,dimnames=list(NULL,profile_columns)))
png(file.path(out,'top-dmrs.png'),width=1400,height=800,res=140)
if (nrow(dmrs)) {
  top <- dmrs[1,,drop=FALSE]
  chr_col <- intersect(c('seqnames','chr','chromosome'),names(top))[1]
  start_col <- intersect(c('start','Start'),names(top))[1]
  end_col <- intersect(c('end','End'),names(top))[1]
  if (any(is.na(c(chr_col,start_col,end_col)))) stop('DMR output lacks chromosome, start, or end')
  same_chr <- sub('^chr','',as.character(probes\$chrom),ignore.case=TRUE) == sub('^chr','',as.character(top[[chr_col]][1]),ignore.case=TRUE)
  in_region <- same_chr & probes\$position >= as.integer(top[[start_col]][1]) & probes\$position <= as.integer(top[[end_col]][1])
  indices <- which(in_region)
  if (length(indices) < 2) stop('top DMR contains fewer than two plotted CpGs')
  positions <- as.integer(probes\$position[indices]); ord <- order(positions); indices <- indices[ord]; positions <- positions[ord]
  grid <- seq(min(positions),max(positions),length.out=max(200,length(unique(positions))))
  colours <- c('#0072B2','#D55E00')
  plotted <- list()
  for (group_index in seq_along(levels(groups))) {
    group_name <- levels(groups)[group_index]
    group_values <- beta[indices,groups == group_name,drop=FALSE]
    means <- rowMeans(group_values)
    standard_error <- apply(group_values,1,sd) / sqrt(ncol(group_values))
    standard_error[!is.finite(standard_error)] <- 0
    lower <- pmax(0,means-1.96*standard_error); upper <- pmin(1,means+1.96*standard_error)
    unique_positions <- sort(unique(positions))
    collapse <- function(values) vapply(unique_positions,function(position) mean(values[positions == position]),numeric(1))
    unique_means <- collapse(means); unique_lower <- collapse(lower); unique_upper <- collapse(upper)
    smooth_values <- function(values) {
      if (length(unique_positions) >= 4) predict(smooth.spline(unique_positions,values,spar=${params.profile_spar}),grid)\$y
      else approx(unique_positions,values,xout=grid,rule=2)\$y
    }
    smooth_mean <- pmin(1,pmax(0,smooth_values(unique_means)))
    smooth_lower <- pmin(smooth_mean,pmax(0,smooth_values(unique_lower)))
    smooth_upper <- pmax(smooth_mean,pmin(1,smooth_values(unique_upper)))
    plotted[[group_index]] <- list(name=group_name,n=ncol(group_values),mean=smooth_mean,lower=smooth_lower,upper=smooth_upper)
    profile <- rbind(profile,data.frame(chrom=as.character(top[[chr_col]][1]),position=grid,group=group_name,mean_beta=approx(unique_positions,unique_means,xout=grid,rule=2)\$y,ci_lower=smooth_lower,ci_upper=smooth_upper,smoothed_beta=smooth_mean))
  }
  plot(grid,plotted[[1]]\$mean,type='n',ylim=c(0,1),xlab=paste0(as.character(top[[chr_col]][1]),' genomic position (GRCh38)'),ylab='Mean methylation fraction',main='Top DMR: normalized group profiles',las=1)
  graphics::grid(col='#E5E5E5',lty=1)
  for (group_index in seq_along(plotted)) {
    item <- plotted[[group_index]]; colour <- colours[group_index]
    polygon(c(grid,rev(grid)),c(item\$lower,rev(item\$upper)),col=adjustcolor(colour,alpha.f=0.18),border=NA)
    lines(grid,item\$mean,col=colour,lwd=3)
  }
  legend('topright',legend=vapply(plotted,function(item) paste0(item\$name,' (n=',item\$n,')'),character(1)),col=colours[seq_along(plotted)],lwd=3,bty='n')
  mtext('Lines are group means smoothed across CpGs; bands are approximate 95% confidence intervals.',side=1,line=4,cex=.8)
} else {
  plot.new(); text(.5,.55,'No DMR passed the configured thresholds',cex=1.2); text(.5,.45,'No regional methylation profile is available',cex=.9,col='#666666')
}
dev.off()
write.table(profile,file.path(out,'top-dmr-profile.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
RS
    """

    stub:
    """
    mkdir dmr
    printf 'probe_id\tchrom\tposition\tdelta_beta\tt\tp_value\tfdr\n' > dmr/cpg-results.tsv
    printf 'seqnames\tstart\tend\n' > dmr/dmr-results.tsv
    printf 'metric\tvalue\ndmrs\t0\n' > dmr/dmr-metrics.tsv
    printf 'chrom\tposition\tgroup\tmean_beta\tci_lower\tci_upper\tsmoothed_beta\n' > dmr/top-dmr-profile.tsv
    touch dmr/top-dmrs.png
    """
}

process COLLECT_RESULTS {
    tag params.cohort_id
    cpus 1
    memory '1 GB'
    container params.python_image

    input:
    path validated
    path dmr

    output:
    path 'results'

    script:
    """
    mkdir results
    cp dmr/cpg-results.tsv results/${params.cohort_id}.cpg-results.tsv
    cp dmr/dmr-results.tsv results/${params.cohort_id}.dmr-results.tsv
    cp dmr/dmr-metrics.tsv results/${params.cohort_id}.dmr-metrics.tsv
    cp dmr/top-dmr-profile.tsv results/${params.cohort_id}.top-dmr-profile.tsv
    cp dmr/top-dmrs.png results/${params.cohort_id}.top-dmrs.png
    cp validated/validation-metrics.tsv results/${params.cohort_id}.validation-metrics.tsv
    printf 'method\tDMRcate kernel smoothing over limma moderated CpG tests\nphenotype\t%s\nreference_genome\tGRCh38\nnot_clinical_use\ttrue\n' '${params.phenotype_column}' > results/${params.cohort_id}.provenance.tsv
    """
}

workflow {
    if (!params.methylation_matrix_uri || !params.phenotype_uri || !params.probe_manifest_uri || !params.phenotype_column || !params.cohort_id) error 'methylation matrix, phenotype, probe manifest, phenotype column, and cohort_id are required'
    if (!(params.cohort_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) error 'invalid cohort_id'
    if ((params.profile_spar as BigDecimal) < 0 || (params.profile_spar as BigDecimal) > 1) error 'profile_spar must be between 0 and 1'
    matrix=channel.fromPath(params.methylation_matrix_uri,checkIfExists:true)
    phenotypes=channel.fromPath(params.phenotype_uri,checkIfExists:true)
    manifest=channel.fromPath(params.probe_manifest_uri,checkIfExists:true)
    validated=VALIDATE_DMR_INPUTS(matrix,phenotypes,manifest)
    dmr=CALL_DMRS(validated)
    results=COLLECT_RESULTS(validated,dmr)
    WRITE_RESULTS(results)
}
