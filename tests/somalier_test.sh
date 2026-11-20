#!/bin/bash
# =============================================================================
# Test of the sex and identity check inside an image: the real somalier on the
# tiny reference and sites built by tests/resources_test.sh. Used by
# tests/ci_checks.sh.
#
#   bash tests/somalier_test.sh IMAGE RESOURCES_DIR
#
# Three samples are written as duplicate-marked BAMs, about 30x at every site:
# A and B with the same genotypes (a duplicate), C with others. The
# somalier_extract and somalier_relate rules run through Snakemake (the BAMs
# and everything upstream are older than them, so nothing else runs), then
# qc_flags.py reads the somalier output. Checked: A and B flagged
# POSSIBLE_DUPLICATE with each other, C not flagged, and the somalier columns
# that qc_flags.py reads are present.
# =============================================================================

set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: bash tests/somalier_test.sh IMAGE RESOURCES_DIR" >&2; exit 1; }
IMAGE="$1"
RESOURCES="$(cd "$2" && pwd)"
IMAGE_HOME="/opt/hemaclone"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "${WORK}/project" "${WORK}/tmp"

docker run --rm -u "$(id -u):$(id -g)" \
    -v "${WORK}:${WORK}" -v "${RESOURCES}:${RESOURCES}" -w "${WORK}/project" -e HOME=/tmp \
    "$IMAGE" bash -c "
set -euo pipefail
mkdir -p bam/units input config qc/markdup
python - <<'PY'
import gzip, random, yaml
random.seed(3)
res = '${RESOURCES}'
ref = {}
name = None
for line in open(res + '/Homo_sapiens_assembly38.fasta'):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:].split()[0]
        ref[name] = []
    else:
        ref[name].append(line)
ref = {k: ''.join(v) for k, v in ref.items()}
sites = []
for line in gzip.open(res + '/sites.hg38.vcf.gz', 'rt'):
    if not line.startswith('#'):
        f = line.split('\t')
        sites.append((int(f[1]), f[3], f[4]))
# genotypes: number of alt alleles (0, 1, 2) per site
shared = [random.choice([0, 1, 1, 2]) for _ in sites]
other = [random.choice([0, 1, 1, 2]) for _ in sites]
genotypes = {'A': shared, 'B': shared, 'C': other}
lengths = {name: len(seq) for name, seq in ref.items()}
for sample, gts in genotypes.items():
    with open(f'{sample}.sam', 'w') as sam:
        sam.write('@HD\tVN:1.6\tSO:unsorted\n')
        for contig, length in lengths.items():
            sam.write(f'@SQ\tSN:{contig}\tLN:{length}\n')
        sam.write(f'@RG\tID:{sample}.B1\tSM:{sample}\tLB:lib\tPL:ILLUMINA\n')
        n = 0
        for (pos, refb, altb), gt in zip(sites, gts):
            for i in range(30):
                start = pos - random.randint(5, 45)
                seq = list(ref['chr1'][start - 1:start - 1 + 50])
                if random.random() < gt / 2:
                    seq[pos - start] = altb
                n += 1
                sam.write(f'r{n}\t0\tchr1\t{start}\t60\t50M\t*\t0\t0\t{\"\".join(seq)}\t' + 'I' * 50 + f'\tRG:Z:{sample}.B1\n')
cfg = yaml.safe_load(open('${IMAGE_HOME}/config/config.yaml'))
settings = {'tools': {'gatk_version': None}, 'reference': {'dir': res}, 'tmpdir': '${WORK}/tmp',
            'job_sizes': {k: {'threads': 1, 'mem_mb': 1000} for k in cfg['job_sizes']}}
yaml.safe_dump(settings, open('${WORK}/settings.yaml', 'w'))
with open('config/samples.tsv', 'w') as sheet:
    sheet.write('sample\tunit\tfq1\tfq2\tflowcell\tlane\tbarcode\trg_id\trg_pu\tlb\tpattern\tbam\n')
    for s in genotypes:
        sheet.write(f'{s}\tB1\t\t\tNA\tNA\tNA\t{s}.B1\t{s}.B1\tlib\tbam\tinput/{s}.bam\n')
PY
# Upstream files first, older, so that only the somalier rules are out of date
for s in A B C; do
    samtools view -b -o input/\$s.bam \$s.sam
    touch -d '2 days ago' input/\$s.bam
    cp input/\$s.bam bam/units/\$s.B1.remapped.bam
    touch -d '1 day ago' bam/units/\$s.B1.remapped.bam
    samtools sort -o bam/\$s.markdup.bam \$s.sam
    samtools index bam/\$s.markdup.bam bam/\$s.markdup.bai
    echo metrics > qc/markdup/\$s.dup_metrics.txt
done
snakemake qc/somalier/relate.samples.tsv qc/somalier/relate.pairs.tsv \
    --snakefile ${IMAGE_HOME}/workflow/Snakefile --directory . \
    --workflow-profile ${IMAGE_HOME}/config/profile \
    --cores 2 --resources mem_mb=3000 \
    --configfile ${WORK}/settings.yaml --config samples=config/samples.tsv \
    --allowed-rules somalier_extract somalier_relate \
    --quiet rules
head -3 qc/somalier/relate.pairs.tsv | cut -f1-4
python - <<'PY'
import sys, types
sys.path.insert(0, '${IMAGE_HOME}/workflow/scripts')
picard = '## METRICS CLASS\tx\nMEAN_COVERAGE\tPCT_20X\tPCT_EXC_DUPE\tPCT_EXC_MAPQ\tMEDIAN_INSERT_SIZE\tPERCENT_DUPLICATION\n30\t0.9\t0.1\t0.01\t300\t0.1\n'
open('metrics.txt', 'w').write(picard)
sm = types.SimpleNamespace(
    input=types.SimpleNamespace(coverage=['metrics.txt'] * 3, insert_size=['metrics.txt'] * 3,
                                dup=['metrics.txt'] * 3,
                                relate_samples='qc/somalier/relate.samples.tsv',
                                relate_pairs='qc/somalier/relate.pairs.tsv'),
    output=['flags.tsv'], log=['flags.log'],
    params=types.SimpleNamespace(samples=['A', 'B', 'C'], data_type='wgs',
                                 thresholds={'min_mean_coverage': 25, 'min_pct_20x': 0.8, 'max_pct_dupe': 0.3},
                                 identity={'min_relatedness_flag': 0.8, 'min_relatedness_list': 0.4,
                                           'sex_y_ratio_male': None, 'expected_sex': None}))
path = '${IMAGE_HOME}/workflow/scripts/qc_flags.py'
try:
    exec(compile(open(path).read(), path, 'exec'), {'snakemake': sm})
except SystemExit as stop:
    sys.exit(f'qc_flags.py stopped: {stop}')
rows = [line.rstrip('\n').split('\t') for line in open('flags.tsv')]
header = rows[0]
table = {r[0]: dict(zip(header, r)) for r in rows[1:]}
for s in 'ABC':
    print(s, table[s]['related_samples'], table[s]['flags'])
assert 'POSSIBLE_DUPLICATE' in table['A']['flags'] and table['A']['related_samples'].startswith('B:'), 'A not flagged with B'
assert 'POSSIBLE_DUPLICATE' in table['B']['flags'], 'B not flagged'
assert 'POSSIBLE_DUPLICATE' not in table['C']['flags'], 'C flagged'
PY
"
echo "Somalier test passed."
