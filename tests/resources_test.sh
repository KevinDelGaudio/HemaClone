#!/bin/bash
# =============================================================================
# Test of 'hemaclone.sh resources' inside an image, on a tiny synthetic bundle
# served from inside the container. Used by tests/ci_checks.sh.
#
#   bash tests/resources_test.sh IMAGE RESOURCES_DIR
#
# Real bwa-mem2, gatk BedToIntervalList, bgzip and tabix run on a 3-contig
# reference with an ALT contig and a GIAB-style repeats BED. The Funcotator step is pre-seeded (the real data sources are tens
# of GB): a folder with a GENCODE-style GTF and its manifest row, so the step
# is skipped as on a re-run. The resulting folder is complete for a run with
# default file names, and the gene panel BED is checked against the expected
# coordinates.
# =============================================================================

set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: bash tests/resources_test.sh IMAGE RESOURCES_DIR" >&2; exit 1; }
IMAGE="$1"
DEST="$2"
[[ ! -e "$DEST" ]] || { echo "ERROR: $DEST already exists" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$DEST"
DEST="$(cd "$DEST" && pwd)"

python3 - "$WORK" "$DEST" <<'EOF'
import hashlib, os, random, sys
work, dest = sys.argv[1:3]
random.seed(1)
bundle = os.path.join(work, "bundle")
somatic = os.path.join(work, "somatic")
os.makedirs(bundle)
os.makedirs(somatic)

contigs = {"chr1": 20000, "chr2": 15000, "chr1_KI270706v1_alt": 500}
with open(os.path.join(bundle, "Homo_sapiens_assembly38.fasta"), "w") as fasta, \
     open(os.path.join(bundle, "Homo_sapiens_assembly38.fasta.fai"), "w") as fai:
    offset = 0
    for name, length in contigs.items():
        fasta.write(f">{name}\n")
        offset += len(name) + 2
        fai.write(f"{name}\t{length}\t{offset}\t60\t61\n")
        seq = "".join(random.choice("ACGT") for _ in range(length))
        for i in range(0, length, 60):
            fasta.write(seq[i:i + 60] + "\n")
        offset += length + (length + 59) // 60
with open(os.path.join(bundle, "Homo_sapiens_assembly38.dict"), "w") as handle:
    handle.write("@HD\tVN:1.6\n")
    for name, length in contigs.items():
        handle.write(f"@SQ\tSN:{name}\tLN:{length}\n")
with open(os.path.join(bundle, "Homo_sapiens_assembly38.fasta.64.alt"), "w") as handle:
    handle.write("chr1_KI270706v1_alt\t0\tchr1\t100\t60\t500M\t*\t0\t0\t*\t*\n")
for name in ["Homo_sapiens_assembly38.dbsnp138.vcf", "Homo_sapiens_assembly38.dbsnp138.vcf.idx",
             "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
             "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
             "wgs_calling_regions.hg38.interval_list"]:
    with open(os.path.join(bundle, name), "wb") as handle:
        handle.write(os.urandom(3000))
for name in ["af-only-gnomad.hg38.vcf.gz", "af-only-gnomad.hg38.vcf.gz.tbi",
             "small_exac_common_3.hg38.vcf.gz", "small_exac_common_3.hg38.vcf.gz.tbi",
             "1000g_pon.hg38.vcf.gz", "1000g_pon.hg38.vcf.gz.tbi"]:
    with open(os.path.join(somatic, name), "wb") as handle:
        handle.write(os.urandom(2000))

# GIAB-style repeats BED, plain gzip as published
import gzip
os.makedirs(os.path.join(work, "giab"))
with gzip.open(os.path.join(work, "giab", "GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz"), "wt") as handle:
    handle.write("chr1\t1995\t2010\nchr1\t5000\t5020\nchr2\t100\t130\n")

# Kit targets: GENEA (two pieces), GENEB, GENED (at the end of chr1); GENEC uncovered
with open(os.path.join(work, "kit.bed"), "w") as handle:
    handle.write("chr1\t950\t1100\tt1\nchr1\t1500\t1600\tt2\nchr2\t4950\t5100\tt3\nchr1\t19900\t20000\tt4\n")

# Pre-seeded Funcotator folder and its manifest row
gencode = os.path.join(dest, "funcotator_dataSources", "gencode", "hg38")
os.makedirs(gencode)
genes = [("chr1", 1001, 1600, "GENEA"), ("chr2", 5001, 5200, "GENEB"),
         ("chr2", 9001, 9500, "GENEC"), ("chr1", 19950, 20000, "GENED")]
with open(os.path.join(gencode, "gencode.v43.annotation.gtf"), "w") as handle:
    for contig, start, end, name in genes:
        handle.write(f'{contig}\tHAVANA\tgene\t{start}\t{end}\t.\t+\t.\t'
                     f'gene_id "{name}.1"; gene_type "protein_coding"; gene_name "{name}";\n')
version = os.path.join(dest, "funcotator_dataSources", "VERSION")
with open(version, "w") as handle:
    handle.write("test-seed\n")
md5 = hashlib.md5(open(version, "rb").read()).hexdigest()
with open(os.path.join(dest, "MANIFEST.tsv"), "w") as handle:
    handle.write("file\tsource\tbytes\tmd5\tdate\n")
    handle.write(f"funcotator_dataSources/VERSION\ttest seed\t{os.path.getsize(version)}\t{md5}\t-\n")
EOF

docker run --rm -u "$(id -u):$(id -g)" \
    -v "${WORK}:${WORK}" -v "${DEST}:${DEST}" -w "$WORK" \
    "$IMAGE" bash -c "
        set -e
        python -m http.server 8765 --bind 127.0.0.1 > server.log 2>&1 &
        sleep 2
        python /opt/hemaclone/workflow/scripts/get_resources.py \
            --resources-dir '$DEST' --targets '${WORK}/kit.bed' \
            --cores 2 --mem-gb 4 --tmpdir '${WORK}/tmp' \
            --set resources.bundle_url=http://127.0.0.1:8765/bundle \
            --set resources.somatic_url=http://127.0.0.1:8765/somatic \
            --set resources.repeats_url=http://127.0.0.1:8765/giab/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz \
            --set 'resources.panel_genes=[GENEA,GENEB,GENEC,GENED]'
    "

EXPECTED=$'chr1\t950\t1150\tGENEA\nchr1\t1450\t1650\tGENEA\nchr1\t19899\t20000\tGENED\nchr2\t4950\t5150\tGENEB'
if [[ "$(cat "${DEST}/chip_genes.bed")" != "$EXPECTED" ]]; then
    echo "ERROR: chip_genes.bed differs from the expected panel:" >&2
    cat "${DEST}/chip_genes.bed" >&2
    exit 1
fi
for name in Homo_sapiens_assembly38.fasta.alt Homo_sapiens_assembly38.fasta.bwt.2bit.64 \
            exome_targets.hg38.interval_list MANIFEST.tsv \
            GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz \
            GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz.tbi; do
    [[ -s "${DEST}/${name}" ]] || { echo "ERROR: ${name} missing or empty" >&2; exit 1; }
done
echo "Resources test passed: $(($(wc -l < "${DEST}/MANIFEST.tsv") - 1)) files in the manifest, panel as expected."
