nextflow.enable.dsl = 2

include { WRITE_RESULTS } from '../base_nf/modules/bundle'

params.genotype_uri = null
params.genotype_format = 'vcf'
params.bgen_sample_uri = ''
params.cohort_id = null
params.dap_output_uri = null
params.reference_dir = '/reference/ancestry/1kgp'
params.reference_release = '1kgp-grch38'
params.n_pcs = 10
params.min_site_match_rate = 0.90
params.min_assignment_probability = 0.80
params.cpus = 8
params.projection_memory = '32 GB'

process PROJECT_ONTO_REFERENCE {
    tag params.cohort_id
    cpus params.cpus
    memory params.projection_memory
    container params.plink2_image

    input:
    path genotype
    path bgen_sample

    output:
    path 'projection'

    script:
    def importArgs = params.genotype_format == 'bgen'
        ? "--bgen '${genotype}' ref-first --sample '${bgen_sample}'"
        : "--vcf '${genotype}'"
    def lastPcColumn = 5 + (params.n_pcs as Integer)
    """
    set -euo pipefail
    test -s '${params.reference_dir}/ref-pc-loadings.tsv.gz'
    test -s '${params.reference_dir}/ref-pc.acount'
    test -s '${params.reference_dir}/ref-pc-scores.tsv.gz'
    mkdir projection
    plink2 ${importArgs} --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 missing --make-pgen --out projection/raw
    gzip -cd '${params.reference_dir}/ref-pc-loadings.tsv.gz' | awk 'NR > 1 {print \$2}' > projection/reference.ids
    plink2 --pfile projection/raw --extract projection/reference.ids --make-pgen --out projection/subset
    reference_sites=\$(wc -l < projection/reference.ids)
    matched_sites=\$(awk 'NR > 1 {n++} END {print n+0}' projection/subset.pvar)
    match_rate=\$(awk -v matched="\$matched_sites" -v total="\$reference_sites" 'BEGIN {if(total<1) exit 2; printf "%.8f", matched/total}')
    if ! awk -v rate="\$match_rate" -v minimum='${params.min_site_match_rate}' 'BEGIN {exit !(rate >= minimum)}'; then
      echo "Reference-site match rate \$match_rate is below ${params.min_site_match_rate}" >&2
      exit 1
    fi
    plink2 --pfile projection/subset --read-freq '${params.reference_dir}/ref-pc.acount' \
      --score '${params.reference_dir}/ref-pc-loadings.tsv.gz' 2 5 header-read no-mean-imputation variance-standardize cols=+scoresums \
      --score-col-nums 6-${lastPcColumn} --out projection/cohort
    printf 'reference_release\t%s\nreference_sites\t%s\nmatched_sites\t%s\nsite_match_rate\t%s\n' \
      '${params.reference_release}' "\$reference_sites" "\$matched_sites" "\$match_rate" > projection/projection-metrics.tsv
    """

    stub:
    """
    mkdir projection
    printf '#FID\tIID\tPC1_AVG\tPC2_AVG\n0\tS1\t0.1\t0.2\n' > projection/cohort.sscore
    printf 'reference_release\tstub\nreference_sites\t100\nmatched_sites\t100\nsite_match_rate\t1\n' > projection/projection-metrics.tsv
    """
}

process CLASSIFY_REFERENCE_SIMILARITY {
    tag params.cohort_id
    cpus 1
    memory '4 GB'
    container params.python_image

    input:
    path projection

    output:
    path 'classification'

    script:
    """
    set -euo pipefail
    mkdir classification
    python3 - projection/cohort.sscore '${params.reference_dir}/ref-pc-scores.tsv.gz' classification ${params.n_pcs} ${params.min_assignment_probability} <<'PY'
import csv, gzip, html, math, pathlib, statistics, sys
projected_path, reference_path, output_dir = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
n_pcs, threshold = int(sys.argv[4]), float(sys.argv[5])
def read_table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
reference, projected = read_table(reference_path), read_table(projected_path)
if not reference or not projected: raise SystemExit('reference and projected PC tables must be non-empty')
label_key = next((key for key in ('population','superpopulation','group') if key in reference[0]), None)
if not label_key: raise SystemExit('reference scores require population, superpopulation, or group')
def pc_keys(row):
    direct=[f'PC{i}' for i in range(1,n_pcs+1)]
    if all(key in row for key in direct): return direct
    averaged=[f'PC{i}_AVG' for i in range(1,n_pcs+1)]
    if all(key in row for key in averaged): return averaged
    scores=[key for key in row if key.endswith('_AVG') and key[:-4].upper().startswith(('PC','SCORE'))]
    return scores[:n_pcs]
ref_keys, sample_keys = pc_keys(reference[0]), pc_keys(projected[0])
if len(ref_keys) < 2 or len(sample_keys) != len(ref_keys): raise SystemExit('reference/projected PC columns do not agree')
ref_vectors=[([float(row[key]) for key in ref_keys], row[label_key]) for row in reference if row.get(label_key)]
means=[statistics.fmean(vector[i] for vector,_ in ref_vectors) for i in range(len(ref_keys))]
sds=[statistics.stdev(vector[i] for vector,_ in ref_vectors) or 1.0 for i in range(len(ref_keys))]
groups=sorted(set(label for _,label in ref_vectors)); centroids={}
for group in groups:
    vectors=[[((value-means[i])/sds[i]) for i,value in enumerate(vector)] for vector,label in ref_vectors if label==group]
    centroids[group]=[statistics.fmean(vector[i] for vector in vectors) for i in range(len(ref_keys))]
rows=[]
for row in projected:
    sample=row.get('IID') or row.get('#IID') or row.get('sample_id') or ''
    vector=[(float(row[key])-means[i])/sds[i] for i,key in enumerate(sample_keys)]
    distances={group:sum((value-centroids[group][i])**2 for i,value in enumerate(vector)) for group in groups}
    weights={group:math.exp(-min(distance,140)/2) for group,distance in distances.items()}; total=sum(weights.values()) or 1
    probabilities={group:value/total for group,value in weights.items()}; best=max(groups,key=probabilities.get); probability=probabilities[best]
    rows.append((sample, best if probability >= threshold else 'unassigned', probability, [float(row[key]) for key in sample_keys]))
with (output_dir/'ancestry.tsv').open('w') as handle:
    handle.write('sample_id\tgenetic_similarity_group\tprobability\t'+'\t'.join(f'PC{i}' for i in range(1,len(ref_keys)+1))+'\n')
    for sample,group,probability,vector in rows: handle.write(f'{sample}\t{group}\t{probability:.6g}\t'+'\t'.join(f'{v:.8g}' for v in vector)+'\n')
counts={group:sum(row[1]==group for row in rows) for group in groups+['unassigned']}
with (output_dir/'classification-metrics.tsv').open('w') as handle:
    handle.write(f'samples\t{len(rows)}\nthreshold\t{threshold}\n'); [handle.write(f'group_{group}\t{count}\n') for group,count in counts.items()]
palette=['#2a78d6','#eb6834','#1baf7a','#8d64d8','#d3a52b','#777']
all_points=[(vector[0],vector[1],label) for vector,label in ref_vectors]+[(vector[0],vector[1],group) for _,group,_,vector in rows]
xs=[p[0] for p in all_points]; ys=[p[1] for p in all_points]; xmin,xmax=min(xs),max(xs); ymin,ymax=min(ys),max(ys)
scale=lambda value,lo,hi,start,size:start+size*(value-lo)/(hi-lo or 1)
svg=[]
for vector,label in ref_vectors: svg.append(f'<circle cx="{scale(vector[0],xmin,xmax,40,720):.1f}" cy="{scale(vector[1],ymax,ymin,20,480):.1f}" r="1.4" fill="#bbb" opacity=".5"/>')
for sample,group,_,vector in rows:
    color=palette[(groups+['unassigned']).index(group)%len(palette)]; svg.append(f'<circle cx="{scale(vector[0],xmin,xmax,40,720):.1f}" cy="{scale(vector[1],ymax,ymin,20,480):.1f}" r="3" fill="{color}"><title>{html.escape(sample)}: {html.escape(group)}</title></circle>')
(output_dir/'pca-overlay.html').write_text('<!doctype html><meta charset="utf-8"><title>1KGP projection</title><h1>Genetic similarity projection</h1><p>Grey: 1KGP reference; coloured: cohort.</p><svg viewBox="0 0 800 520" width="100%">'+''.join(svg)+'</svg>')
PY
    """

    stub:
    """
    mkdir classification
    printf 'sample_id\tgenetic_similarity_group\tprobability\tPC1\tPC2\nS1\tEAS\t0.99\t0.1\t0.2\n' > classification/ancestry.tsv
    printf 'samples\t1\nthreshold\t0.8\ngroup_EAS\t1\n' > classification/classification-metrics.tsv
    touch classification/pca-overlay.html
    """
}

process COLLECT_RESULTS {
    tag params.cohort_id
    cpus 1
    memory '1 GB'
    container params.python_image

    input:
    path projection
    path classification

    output:
    path 'results'

    script:
    """
    mkdir results
    cp classification/ancestry.tsv results/${params.cohort_id}.ancestry.tsv
    cp classification/pca-overlay.html results/${params.cohort_id}.pca-overlay.html
    cp classification/classification-metrics.tsv results/${params.cohort_id}.classification-metrics.tsv
    cp projection/projection-metrics.tsv results/${params.cohort_id}.projection-metrics.tsv
    printf 'reference_release\t%s\nreference_dir\t%s\ninterpretation\tgenetic similarity relative to 1KGP; not ethnicity or nationality\n' '${params.reference_release}' '${params.reference_dir}' > results/${params.cohort_id}.provenance.tsv
    """
}

workflow {
    if (!params.genotype_uri || !params.cohort_id) error 'genotype_uri and cohort_id are required'
    if (!(params.genotype_format in ['vcf','bgen'])) error 'genotype_format must be vcf or bgen'
    if (params.genotype_format == 'bgen' && !params.bgen_sample_uri) error 'bgen_sample_uri is required for BGEN input'
    if (!(params.cohort_id ==~ /[A-Za-z0-9][A-Za-z0-9._-]{0,127}/)) error 'invalid cohort_id'
    genotype=channel.fromPath(params.genotype_uri, checkIfExists:true)
    bgen_sample=params.bgen_sample_uri ? channel.fromPath(params.bgen_sample_uri, checkIfExists:true) : channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists:true)
    projection=PROJECT_ONTO_REFERENCE(genotype,bgen_sample)
    classification=CLASSIFY_REFERENCE_SIMILARITY(projection)
    results=COLLECT_RESULTS(projection,classification)
    WRITE_RESULTS(results)
}
