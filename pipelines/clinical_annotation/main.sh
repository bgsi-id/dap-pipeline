#!/usr/bin/env bash
#
# annotate_wgs.sh - fast WGS annotation for rare disease
#
# Flow: normalize -> echtvar AF/ClinVar join -> prefilter -> VEP on reduced set
#
# Usage:
#   ./annotate_wgs.sh -i sample.vcf.gz -s SAMPLE01 -o out/ -r /home/ubuntu/reference
#
set -euo pipefail

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
THREADS=16
AF_CUTOFF=0.01
FASTA_NAME="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
VEP_CACHE_SUBDIR="vep"

IMG_BCFTOOLS="quay.io/biocontainers/bcftools:1.21--h8b25389_0"
IMG_ECHTVAR="quay.io/biocontainers/echtvar:0.2.2--h4349ce8_0"
IMG_VEP="ensemblorg/ensembl-vep:release_116.0"

# ----------------------------------------------------------------------------
# args
# ----------------------------------------------------------------------------
IN=""; SAMPLE=""; OUTDIR=""; REFDIR=""
while getopts "i:s:o:r:t:f:h" opt; do
  case $opt in
    i) IN=$(readlink -f "$OPTARG") ;;
    s) SAMPLE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    r) REFDIR=$(readlink -f "$OPTARG") ;;
    t) THREADS="$OPTARG" ;;
    f) AF_CUTOFF="$OPTARG" ;;
    h) sed -n '2,10p' "$0"; exit 0 ;;
    *) exit 1 ;;
  esac
done

[[ -z "$IN" || -z "$SAMPLE" || -z "$OUTDIR" || -z "$REFDIR" ]] && {
  echo "usage: $0 -i <vcf.gz> -s <sample> -o <outdir> -r <refdir> [-t threads] [-f af_cutoff]" >&2
  exit 1
}

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")
WORK="$OUTDIR/work"
mkdir -p "$WORK"
LOG="$OUTDIR/${SAMPLE}.log"

UID_GID="$(id -u):$(id -g)"
IN_DIR=$(dirname "$IN")
IN_BASE=$(basename "$IN")

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# docker wrappers. /ref, /work, /in are the mount points inside every container
drun() {
  local img="$1"; shift
  docker run --rm --user "$UID_GID" \
    -v "$REFDIR":/ref:ro \
    -v "$WORK":/work \
    -v "$IN_DIR":/in:ro \
    -v "$OUTDIR":/out \
    -w /work \
    "$img" "$@"
}

bcf()     { drun "$IMG_BCFTOOLS" "$@"; }
echtvar() { drun "$IMG_ECHTVAR" echtvar "$@"; }
vep()     { drun "$IMG_VEP" "$@"; }

# ----------------------------------------------------------------------------
# preflight
# ----------------------------------------------------------------------------
log "preflight"
for f in "$FASTA_NAME" "gnomad_v4.1.zip" "clinvar.chr.vcf.gz" "$VEP_CACHE_SUBDIR"; do
  [[ -e "$REFDIR/$f" ]] || { echo "MISSING: $REFDIR/$f" >&2; exit 1; }
done
command -v docker >/dev/null || { echo "docker not found" >&2; exit 1; }

SPLICEAI="spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz"
USE_SPLICEAI=0
if [[ -s "$REFDIR/$SPLICEAI" && -s "$REFDIR/$SPLICEAI.tbi" ]]; then
  USE_SPLICEAI=1
else
  log "WARN: SpliceAI file absent or unindexed, skipping splice annotation"
fi

# ----------------------------------------------------------------------------
# 1. normalize
# ----------------------------------------------------------------------------
log "1/5 normalize"
bcf bash -c "
  bcftools norm -m -any -f /ref/$FASTA_NAME -Ou /in/$IN_BASE \
  | bcftools norm -d exact --threads $THREADS -Oz -o /work/norm.vcf.gz
  bcftools index -t --threads $THREADS /work/norm.vcf.gz
" 2>&1 | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 2. frequency and ClinVar join
# ----------------------------------------------------------------------------
log "2/5 echtvar gnomAD join"
echtvar anno -e /ref/gnomad_v4.1.zip /work/norm.vcf.gz /work/freq.vcf.gz 2>&1 | tee -a "$LOG"

log "2b/5 ClinVar join"
bcf bash -c "
  bcftools index -t --threads $THREADS /work/freq.vcf.gz
  bcftools annotate -a /ref/clinvar.chr.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVI \
    --threads $THREADS -Oz -o /work/annot.vcf.gz /work/freq.vcf.gz
  bcftools index -t --threads $THREADS /work/annot.vcf.gz
" 2>&1 | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 3. prefilter. rare OR ClinVar pathogenic
# ----------------------------------------------------------------------------
log "3/5 prefilter (AF < $AF_CUTOFF OR ClinVar P/LP)"
bcf bash -c "
  bcftools filter -i '
    ((INFO/gnomad_af < $AF_CUTOFF || INFO/gnomad_af = -1) &&
     (INFO/gnomad_af_max < 0.02 || INFO/gnomad_af_max = -1))
    || INFO/CLNSIG ~ \"athogenic\"
  ' --threads $THREADS -Oz -o /work/tiered.vcf.gz /work/annot.vcf.gz
  bcftools index -t --threads $THREADS /work/tiered.vcf.gz
" 2>&1 | tee -a "$LOG"

N_IN=$(bcf bcftools index -n /work/annot.vcf.gz | tr -d '\r')
N_OUT=$(bcf bcftools index -n /work/tiered.vcf.gz | tr -d '\r')
log "    $N_IN -> $N_OUT variants"

# ----------------------------------------------------------------------------
# 4. VEP on the reduced set
# ----------------------------------------------------------------------------
log "4/5 VEP"
VEP_ARGS=(
  vep
  --input_file /work/tiered.vcf.gz
  --output_file /out/${SAMPLE}.annotated.vcf.gz
  --vcf --compress_output bgzip
  --offline --cache --dir_cache /ref/$VEP_CACHE_SUBDIR --merged
  --cache_version 116 --assembly GRCh38
  --fasta /ref/$FASTA_NAME
  --mane --mane_select --canonical --symbol --biotype
  --hgvs --hgvsg --shift_hgvs 1
  --numbers --domains --protein --uniprot
  --pick_order mane_select,mane_plus_clinical,canonical,rank
  --fork "$THREADS" --buffer_size 50000
  --no_stats --force_overwrite
)
[[ $USE_SPLICEAI -eq 1 ]] && VEP_ARGS+=( --plugin "SpliceAI,snv=/ref/$SPLICEAI" )

vep "${VEP_ARGS[@]}" 2>&1 | tee -a "$LOG"

# ----------------------------------------------------------------------------
# 5. flat TSV for review
# ----------------------------------------------------------------------------
log "5/5 TSV export"
bcf bash -c "
  bcftools index -t --threads $THREADS /out/${SAMPLE}.annotated.vcf.gz
  bcftools +split-vep -d -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%Feature\t%Consequence\t%IMPACT\t%HGVSc\t%HGVSp\t%MANE_SELECT\t%INFO/gnomad_af\t%INFO/gnomad_af_max\t%INFO/CLNSIG\t%INFO/CLNDN\t[%GT\t%DP\t%GQ]\n' \
    -i 'IMPACT=\"HIGH\" || IMPACT=\"MODERATE\" || CLNSIG ~ \"athogenic\"' \
    /out/${SAMPLE}.annotated.vcf.gz > /out/${SAMPLE}.report.tsv
" 2>&1 | tee -a "$LOG"

sed -i '1i CHROM\tPOS\tREF\tALT\tSYMBOL\tTRANSCRIPT\tCONSEQUENCE\tIMPACT\tHGVSc\tHGVSp\tMANE\tGNOMAD_AF\tGNOMAD_AF_MAX\tCLNSIG\tCLNDN\tGT\tDP\tGQ' \
  "$OUTDIR/${SAMPLE}.report.tsv"

# ----------------------------------------------------------------------------
# provenance
# ----------------------------------------------------------------------------
cat > "$OUTDIR/${SAMPLE}.provenance.txt" <<EOF
sample            $SAMPLE
input             $IN
date              $(date -Iseconds)
af_cutoff         $AF_CUTOFF
reference         $FASTA_NAME
vep               release 116, merged cache
mane              v1.5 summary (report layer)
gnomad            v4.1 genomes, AF and AF_grpmax only
clinvar           $(ls -l --time-style=+%Y-%m-%d "$REFDIR/clinvar.chr.vcf.gz" | awk '{print $6}')
spliceai          $([[ $USE_SPLICEAI -eq 1 ]] && echo "ensembl_mane_v1.4" || echo "NOT APPLIED")
images            $IMG_BCFTOOLS
                  $IMG_ECHTVAR
                  $IMG_VEP
variants_in       $N_IN
variants_tiered   $N_OUT

NOT VALIDATED FOR CLINICAL USE.
Missing: nhomalt (no recessive homozygote filter), internal AF panel,
SV/CNV/repeat annotation, constraint metrics, GIAB concordance run.
EOF

log "done"
log "  $OUTDIR/${SAMPLE}.annotated.vcf.gz"
log "  $OUTDIR/${SAMPLE}.report.tsv"
log "  $OUTDIR/${SAMPLE}.provenance.txt"