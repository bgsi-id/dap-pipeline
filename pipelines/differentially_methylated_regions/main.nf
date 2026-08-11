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
png(file.path(out,'top-dmrs.png'),width=1400,height=800,res=140)
if (nrow(dmrs)) {
  top <- head(dmrs,20); score <- -log10(pmax(top\$Stouffer,1e-300)); barplot(rev(score),names.arg=rev(paste0(top\$seqnames,':',top\$start,'-',top\$end)),horiz=TRUE,las=1,xlab='-log10(Stouffer p)',main='Top differentially methylated regions')
} else plot.new()
dev.off()
RS
    """

    stub:
    """
    mkdir dmr
    printf 'probe_id\tchrom\tposition\tdelta_beta\tt\tp_value\tfdr\n' > dmr/cpg-results.tsv
    printf 'seqnames\tstart\tend\n' > dmr/dmr-results.tsv
    printf 'metric\tvalue\ndmrs\t0\n' > dmr/dmr-metrics.tsv
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
    cp dmr/top-dmrs.png results/${params.cohort_id}.top-dmrs.png
    cp validated/validation-metrics.tsv results/${params.cohort_id}.validation-metrics.tsv
    printf 'method\tDMRcate kernel smoothing over limma moderated CpG tests\nphenotype\t%s\nreference_genome\tGRCh38\nnot_clinical_use\ttrue\n' '${params.phenotype_column}' > results/${params.cohort_id}.provenance.tsv
    """
}

workflow {
    if (!params.methylation_matrix_uri || !params.phenotype_uri || !params.probe_manifest_uri || !params.phenotype_column || !params.cohort_id) error 'methylation matrix, phenotype, probe manifest, phenotype column, and cohort_id are required'
    if (!(params.cohort_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) error 'invalid cohort_id'
    matrix=channel.fromPath(params.methylation_matrix_uri,checkIfExists:true)
    phenotypes=channel.fromPath(params.phenotype_uri,checkIfExists:true)
    manifest=channel.fromPath(params.probe_manifest_uri,checkIfExists:true)
    validated=VALIDATE_DMR_INPUTS(matrix,phenotypes,manifest)
    dmr=CALL_DMRS(validated)
    results=COLLECT_RESULTS(validated,dmr)
    WRITE_RESULTS(results)
}
