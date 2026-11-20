#!/bin/bash
# =============================================================================
# Test of BAM and CRAM input inside an image: real samtools and bwa-mem2 on the
# tiny reference built by tests/resources_test.sh. Used by tests/ci_checks.sh.
#
#   bash tests/bam_input_test.sh IMAGE RESOURCES_DIR
#
# Read pairs simulated from the reference are written as an unaligned BAM and
# as a CRAM (with original qualities in the OQ tag, and one unpaired read);
# the sample sheet is built from them, both files pass bam_check, and
# bam_remap re-aligns them. Checked: every pair re-aligned and mapped, the
# unpaired read dropped, coordinate-sorted output with the expected read group
# (library taken from the header), and the qualities taken from OQ.
# =============================================================================

set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: bash tests/bam_input_test.sh IMAGE RESOURCES_DIR" >&2; exit 1; }
IMAGE="$1"
RESOURCES="$(cd "$2" && pwd)"
IMAGE_HOME="/opt/hemaclone"
PAIRS=50

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "${WORK}/project/fastq" "${WORK}/tmp"

cat > "${WORK}/settings.yaml" <<EOF
tools:
  gatk_version: null
input:
  cram_reference: ${RESOURCES}/Homo_sapiens_assembly38.fasta
reference:
  dir: ${RESOURCES}
tmpdir: ${WORK}/tmp
EOF
docker run --rm -u "$(id -u):$(id -g)" \
    -v "${WORK}:${WORK}" -v "${RESOURCES}:${RESOURCES}" -w "${WORK}/project" -e HOME=/tmp \
    "$IMAGE" bash -c "
set -euo pipefail
python - <<'PY'
import random, yaml
random.seed(2)
ref = {}
name = None
for line in open('${RESOURCES}/Homo_sapiens_assembly38.fasta'):
    line = line.strip()
    if line.startswith('>'):
        name = line[1:].split()[0]
        ref[name] = []
    else:
        ref[name].append(line)
chr1 = ''.join(ref['chr1'])
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(s):
    return ''.join(comp[b] for b in reversed(s))
with open('reads.sam', 'w') as sam:
    sam.write('@HD\tVN:1.6\tSO:unsorted\n@RG\tID:orig\tSM:different_name\tLB:libX\tPL:ILLUMINA\n')
    for i in range($PAIRS):
        start = random.randint(1000, 18000)
        r1 = chr1[start:start + 100]
        r2 = revcomp(chr1[start + 200:start + 300])
        # the read qualities are low (#); the original ones in OQ are high (I)
        tags = 'RG:Z:orig\tOQ:Z:' + 'I' * 100
        sam.write(f'p{i}\t77\t*\t0\t0\t*\t*\t0\t0\t{r1}\t' + '#' * 100 + f'\t{tags}\n')
        sam.write(f'p{i}\t141\t*\t0\t0\t*\t*\t0\t0\t{r2}\t' + '#' * 100 + f'\t{tags}\n')
    sam.write('lonely\t4\t*\t0\t0\t*\t*\t0\t0\t' + chr1[5000:5100] + '\t' + 'I' * 100 + '\tRG:Z:orig\n')
s = yaml.safe_load(open('${WORK}/settings.yaml'))
cfg = yaml.safe_load(open('${IMAGE_HOME}/config/config.yaml'))
s['job_sizes'] = {k: {'threads': 1, 'mem_mb': 1000} for k in cfg['job_sizes']}
yaml.safe_dump(s, open('${WORK}/settings.yaml', 'w'))
PY
samtools view -b -o fastq/SIMD.bam reads.sam
samtools view -C -T '${RESOURCES}/Homo_sapiens_assembly38.fasta' -o fastq/SIME.cram reads.sam
python ${IMAGE_HOME}/workflow/scripts/make_samplesheet.py --fastq-dir fastq --out config/samples.tsv
snakemake qc/bam_check/SIMD.B1.ok qc/bam_check/SIME.B1.ok \
    bam/units/SIMD.B1.remapped.bam bam/units/SIME.B1.remapped.bam \
    --snakefile ${IMAGE_HOME}/workflow/Snakefile --directory . \
    --workflow-profile ${IMAGE_HOME}/config/profile \
    --cores 2 --resources mem_mb=3000 \
    --configfile ${WORK}/settings.yaml --config samples=config/samples.tsv \
    --notemp --quiet rules
for sample in SIMD SIME; do
    out=bam/units/\${sample}.B1.remapped.bam
    total=\$(samtools view -c \$out)
    mapped=\$(samtools view -c -F 4 \$out)
    sorted=\$(samtools view -H \$out | grep -c 'SO:coordinate' || true)
    rg=\$(samtools view -H \$out | grep '^@RG' || true)
    high=\$(samtools view \$out | awk '\$11 ~ /^I+\$/' | wc -l)
    echo \"\$sample: \$total reads, \$mapped mapped, sorted=\$sorted, reads with OQ qualities=\$high\"
    [ \"\$total\" -eq $((PAIRS * 2)) ] || { echo \"ERROR: \$sample: expected $((PAIRS * 2)) reads (unpaired read dropped)\" >&2; exit 1; }
    [ \"\$mapped\" -eq \"\$total\" ] || { echo \"ERROR: \$sample: unmapped reads\" >&2; exit 1; }
    [ \"\$sorted\" -eq 1 ] || { echo \"ERROR: \$sample: not coordinate-sorted\" >&2; exit 1; }
    [ \"\$high\" -eq \"\$total\" ] || { echo \"ERROR: \$sample: qualities not taken from OQ\" >&2; exit 1; }
    echo \"\$rg\" | grep -q \"ID:\${sample}.B1\" && echo \"\$rg\" | grep -q 'SM:'\${sample} && echo \"\$rg\" | grep -q 'LB:libX' \
        || { echo \"ERROR: \$sample: read group \$rg\" >&2; exit 1; }
done
grep -h 'singletons\|processed' logs/bam_remap/*.log || true
"
echo "BAM/CRAM input test passed."
