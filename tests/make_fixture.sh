#!/bin/bash
# =============================================================================
# Build a minimal synthetic project for dry runs and linting (used by CI).
#
#   bash tests/make_fixture.sh PROJECT_DIR RESOURCES_DIR
#
# PROJECT_DIR gets fastq/ with three tiny gzipped FASTQ pairs: one sample split
# over two lanes and one sample in a single pair, so both layouts are planned.
# RESOURCES_DIR gets empty files under the default resource names, enough for
# the preflight checks, which test that files exist. Nothing here is real data:
# the fixture is for planning the jobs, not for running them.
# =============================================================================

set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: bash tests/make_fixture.sh PROJECT_DIR RESOURCES_DIR" >&2; exit 1; }
PROJECT_DIR="$1"
RESOURCES_DIR="$2"

[[ ! -e "${PROJECT_DIR}/fastq" ]] || { echo "ERROR: ${PROJECT_DIR}/fastq already exists" >&2; exit 1; }
[[ ! -e "$RESOURCES_DIR" ]] || { echo "ERROR: ${RESOURCES_DIR} already exists" >&2; exit 1; }

# One FASTQ file: write_fastq PATH FLOWCELL LANE MATE
write_fastq() {
    local path="$1" flowcell="$2" lane="$3" mate="$4" i
    for i in 1 2 3 4; do
        printf '@SIM01:1:%s:%s:1101:%d:1000 %s:N:0:ACGTACGT+TGCATGCA\n' "$flowcell" "$lane" "$i" "$mate"
        printf 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\n'
        printf 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n'
    done | gzip -c > "$path"
}

mkdir -p "${PROJECT_DIR}/fastq"
for mate in 1 2; do
    write_fastq "${PROJECT_DIR}/fastq/SIMA_S1_L001_R${mate}_001.fastq.gz" FCSIM0001 1 "$mate"
    write_fastq "${PROJECT_DIR}/fastq/SIMA_S1_L002_R${mate}_001.fastq.gz" FCSIM0001 2 "$mate"
    write_fastq "${PROJECT_DIR}/fastq/SIMB_R${mate}.fastq.gz" FCSIM0002 1 "$mate"
done

mkdir -p "${RESOURCES_DIR}/funcotator_dataSources"
for name in \
    Homo_sapiens_assembly38.fasta \
    Homo_sapiens_assembly38.fasta.fai \
    Homo_sapiens_assembly38.dict \
    Homo_sapiens_assembly38.fasta.bwt.2bit.64 \
    exome_targets.hg38.interval_list \
    wgs_calling_regions.hg38.interval_list \
    chip_genes.bed \
    af-only-gnomad.hg38.vcf.gz af-only-gnomad.hg38.vcf.gz.tbi \
    small_exac_common_3.hg38.vcf.gz small_exac_common_3.hg38.vcf.gz.tbi \
    Homo_sapiens_assembly38.dbsnp138.vcf Homo_sapiens_assembly38.dbsnp138.vcf.idx \
    Mills_and_1000G_gold_standard.indels.hg38.vcf.gz Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
    1000g_pon.hg38.vcf.gz 1000g_pon.hg38.vcf.gz.tbi \
    GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz \
    GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz.tbi \
    sites.hg38.vcf.gz
do
    : > "${RESOURCES_DIR}/${name}"
done

echo "FASTQ files : $(ls "${PROJECT_DIR}/fastq" | wc -l) in ${PROJECT_DIR}/fastq"
echo "Resources   : $(ls "$RESOURCES_DIR" | wc -l) entries in ${RESOURCES_DIR}"
