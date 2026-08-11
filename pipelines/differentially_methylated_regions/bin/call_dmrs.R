#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 10) stop('expected matrix_dir, output_dir, partition and seven analysis parameters')
matrix_dir <- args[[1]]
out <- args[[2]]
partition <- args[[3]]
fdr_threshold <- as.numeric(args[[4]])
lambda <- as.numeric(args[[5]])
bandwidth_scaling <- as.numeric(args[[6]])
min_cpgs <- as.integer(args[[7]])
min_delta_beta <- as.numeric(args[[8]])
profile_spar <- as.numeric(args[[9]])
min_sites <- as.integer(args[[10]])

required <- c('limma','DMRcate','GenomicRanges','IRanges')
missing <- required[!vapply(required,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))]
if (length(missing)) stop('DMR image is missing packages: ',paste(missing,collapse=', '))

dir.create(out)
beta_frame <- read.delim(gzfile(file.path(matrix_dir,'beta.tsv.gz')),check.names=FALSE)
probes <- read.delim(gzfile(file.path(matrix_dir,'probes.tsv.gz')),check.names=FALSE)
metadata <- read.delim(file.path(matrix_dir,'metadata.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
if (nrow(beta_frame) < min_sites) stop('too few coverage-filtered CpGs for DMR analysis')
if (!identical(beta_frame$probe_id,probes$probe_id)) stop('beta matrix and probe coordinates disagree')
if (!identical(names(beta_frame)[-1],metadata$sample_id)) stop('beta matrix and cohort membership disagree')
groups <- factor(metadata$group,levels=c('control','patient'))
if (anyNA(groups) || any(table(groups) < 3)) stop('control and patient each require at least three samples')
beta <- as.matrix(beta_frame[,-1,drop=FALSE]); storage.mode(beta) <- 'double'; rownames(beta) <- beta_frame$probe_id
if (any(!is.finite(beta),na.rm=TRUE) || any(beta < 0 | beta > 1,na.rm=TRUE)) stop('methylation fractions must be in [0,1]')

design <- model.matrix(~groups)
clipped <- pmin(pmax(beta,1e-6),1-1e-6)
mvalues <- log2(clipped/(1-clipped))
fit <- limma::eBayes(limma::lmFit(mvalues,design),robust=TRUE)
p <- fit$p.value[,2]
q <- p.adjust(p,'BH')
statistic <- fit$t[,2]
delta <- rowMeans(beta[,groups=='patient',drop=FALSE],na.rm=TRUE)-rowMeans(beta[,groups=='control',drop=FALSE],na.rm=TRUE)
ranges <- GenomicRanges::GRanges(
  seqnames=probes$chrom,
  ranges=IRanges::IRanges(as.integer(probes$position),as.integer(probes$position)),
  stat=statistic,rawpval=p,diff=delta,ind.fdr=q,is.sig=q<fdr_threshold
)
names(ranges) <- rownames(beta)
annotated <- methods::new('CpGannotated',ranges=ranges)
called <- DMRcate::dmrcate(
  annotated,lambda=lambda,C=bandwidth_scaling,min.cpgs=min_cpgs,
  betacutoff=min_delta_beta
)
dmrs <- as.data.frame(DMRcate::extractRanges(called,genome='hg38'))
write.table(
  data.frame(probe_id=rownames(beta),chrom=probes$chrom,position=probes$position,
             delta_beta=delta,t=statistic,p_value=p,fdr=q),
  file.path(out,'cpg-results.tsv'),sep='\t',quote=FALSE,row.names=FALSE
)
write.table(dmrs,file.path(out,'dmr-results.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
metrics <- data.frame(
  metric=c('partition','samples','control_samples','patient_samples','cpgs_tested','cpgs_fdr','dmrs'),
  value=c(partition,ncol(beta),sum(groups=='control'),sum(groups=='patient'),nrow(beta),sum(q<fdr_threshold,na.rm=TRUE),nrow(dmrs))
)
write.table(metrics,file.path(out,'dmr-metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

profile_columns <- c('partition','chrom','position','group','mean_beta','ci_lower','ci_upper','smoothed_beta')
profile <- data.frame(matrix(ncol=length(profile_columns),nrow=0,dimnames=list(NULL,profile_columns)))
png(file.path(out,'top-dmr-profile.png'),width=1400,height=800,res=140)
if (nrow(dmrs)) {
  top <- dmrs[1,,drop=FALSE]
  chr_col <- intersect(c('seqnames','chr','chromosome'),names(top))[1]
  start_col <- intersect(c('start','Start'),names(top))[1]
  end_col <- intersect(c('end','End'),names(top))[1]
  if (any(is.na(c(chr_col,start_col,end_col)))) stop('DMR output lacks chromosome, start, or end')
  same_chr <- sub('^chr','',as.character(probes$chrom),ignore.case=TRUE) == sub('^chr','',as.character(top[[chr_col]][1]),ignore.case=TRUE)
  in_region <- same_chr & probes$position >= as.integer(top[[start_col]][1]) & probes$position <= as.integer(top[[end_col]][1])
  indices <- which(in_region)
  if (length(indices) < 2) stop('top DMR contains fewer than two plotted CpGs')
  positions <- as.integer(probes$position[indices])
  ord <- order(positions); indices <- indices[ord]; positions <- positions[ord]
  plot_x <- seq(min(positions),max(positions),length.out=max(200,length(unique(positions))))
  colours <- c('#0072B2','#D55E00')
  plotted <- list()
  for (group_index in seq_along(levels(groups))) {
    group_name <- levels(groups)[group_index]
    group_values <- beta[indices,groups == group_name,drop=FALSE]
    observations <- rowSums(!is.na(group_values))
    means <- rowMeans(group_values,na.rm=TRUE)
    standard_error <- apply(group_values,1,sd,na.rm=TRUE) / sqrt(observations)
    standard_error[!is.finite(standard_error)] <- 0
    lower <- pmax(0,means-1.96*standard_error)
    upper <- pmin(1,means+1.96*standard_error)
    unique_positions <- sort(unique(positions))
    collapse <- function(values) vapply(unique_positions,function(position) mean(values[positions == position],na.rm=TRUE),numeric(1))
    unique_means <- collapse(means); unique_lower <- collapse(lower); unique_upper <- collapse(upper)
    smooth_values <- function(values) {
      if (length(unique_positions) >= 4) predict(smooth.spline(unique_positions,values,spar=profile_spar),plot_x)$y
      else approx(unique_positions,values,xout=plot_x,rule=2)$y
    }
    smooth_mean <- pmin(1,pmax(0,smooth_values(unique_means)))
    smooth_lower <- pmin(smooth_mean,pmax(0,smooth_values(unique_lower)))
    smooth_upper <- pmax(smooth_mean,pmin(1,smooth_values(unique_upper)))
    plotted[[group_index]] <- list(name=group_name,n=ncol(group_values),mean=smooth_mean,lower=smooth_lower,upper=smooth_upper)
    profile <- rbind(profile,data.frame(
      partition=partition,chrom=as.character(top[[chr_col]][1]),position=plot_x,
      group=group_name,mean_beta=approx(unique_positions,unique_means,xout=plot_x,rule=2)$y,
      ci_lower=smooth_lower,ci_upper=smooth_upper,smoothed_beta=smooth_mean
    ))
  }
  plot(
    plot_x,plotted[[1]]$mean,type='n',ylim=c(0,1),
    xlab=paste0(as.character(top[[chr_col]][1]),' genomic position (GRCh38)'),
    ylab='Mean methylation fraction',
    main=paste0('Top DMR — ',partition,' normalized group profiles'),las=1
  )
  graphics::grid(col='#E5E5E5',lty=1)
  for (group_index in seq_along(plotted)) {
    item <- plotted[[group_index]]; colour <- colours[group_index]
    polygon(c(plot_x,rev(plot_x)),c(item$lower,rev(item$upper)),col=adjustcolor(colour,alpha.f=0.18),border=NA)
    lines(plot_x,item$mean,col=colour,lwd=3)
  }
  legend('topright',legend=vapply(plotted,function(item) paste0(item$name,' (n=',item$n,')'),character(1)),col=colours,lwd=3,bty='n')
  mtext('Group means smoothed across CpGs; bands are approximate 95% confidence intervals.',side=1,line=4,cex=.8)
} else {
  plot.new(); text(.5,.55,paste('No',partition,'DMR passed the configured thresholds'),cex=1.2)
  text(.5,.45,'No regional methylation profile is available',cex=.9,col='#666666')
}
dev.off()
write.table(profile,file.path(out,'top-dmr-profile.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

writeLines(c(
  paste('partition',partition,sep='\t'),
  'contrast\tpatient - control',
  'method\tDMRcate over limma moderated CpG tests',
  'input\tmodkit bedMethyl valid coverage and modified counts',
  'reference_genome\tGRCh38',
  if (partition == 'combined') 'interpretation\tprimary all-read result' else 'interpretation\texploratory haplotype/unphased partition; HP labels are not aligned across samples',
  'not_clinical_use\ttrue'
),file.path(out,'provenance.tsv'))
