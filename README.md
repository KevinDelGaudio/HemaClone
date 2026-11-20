# HemaClone

HemaClone calls clonal haematopoiesis (CH) variants from blood exome
sequencing without a matched normal, and classifies each sample as CHIP+ or
CHIP-. It runs end to end, from FASTQ files to maftools-ready MAFs, in a
single versioned Docker image with hard CPU and memory limits.

- Alignment (bwa-mem2), duplicate marking and base quality recalibration (GATK4)
- Tumor-only Mutect2 on a CH gene panel, with read orientation and contamination models
- External panel of normals (flag or filter), with known CH hotspots exempt
- Allele-matched gnomAD frequencies, Funcotator annotation
- CH classification with one pass/fail column per criterion: protein-altering
  class, alternate reads, depth, VAF window, binomial test against the
  heterozygous fraction, gnomAD AF, panel of normals
- Cohort tables and MAFs for maftools, CHIP- samples included
- QC: FastQC, alignment, insert size and hybrid-selection metrics, MultiQC, per-sample flags

Design, criteria and outputs in detail: [`docs/CONTEXT.md`](docs/CONTEXT.md).
Changes between versions: [`CHANGELOG.md`](CHANGELOG.md). How to contribute: [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Requirements

- Linux with Docker; the user in the `docker` group
- hg38 resources in one folder, built with `bash hemaclone.sh resources` (below)
- The target BED of the exome kit, hg38 coordinates with `chr` names

## Get the image

Pull a released image and give it the local name the launcher expects:

```bash
docker pull ghcr.io/kevindelgaudio/hemaclone:1.6.0
docker tag ghcr.io/kevindelgaudio/hemaclone:1.6.0 hemaclone:1.6.0
```

Or build it from a checkout:

```bash
git clone https://github.com/KevinDelGaudio/HemaClone.git
cd HemaClone
bash workflow/hemaclone.sh build          # builds hemaclone:<VERSION>
```

## Resources

One command downloads and builds everything a run needs, with the file names
the defaults expect:

```bash
docker run --rm hemaclone:1.6.0 cat /opt/hemaclone/workflow/hemaclone.sh > hemaclone.sh

bash hemaclone.sh resources \
    --resources-dir /path/to/hg38_resources \
    --targets /path/to/exome_kit_targets.hg38.bed \
    --cores 8 --mem-gb 96 \
    --tmpdir /path/to/tmp
```

- GATK hg38 reference (with its ALT contig list for alt-aware mapping), dbSNP,
  Mills indels, gnomAD AF-only, common sites and the 1000 Genomes panel of
  normals, from GATK's public buckets; downloads resume and are checked
  against the server's size and MD5.
- GIAB tandem-repeat and homopolymer BED (Dwarshuis et al. 2024), for the
  stricter thresholds on indels in repeats.
- bwa-mem2 index (needs much memory: if it is killed, raise `--mem-gb`),
  Funcotator somatic data sources.
- The kit BED as an interval list, and the CH gene panel BED: the genes in
  `resources.panel_genes`, from the Funcotator GENCODE annotation, intersected
  with the kit targets and padded by `resources.panel_padding` (50 bp). Genes
  the kit does not cover are listed.
- `MANIFEST.tsv` records every file with its source, size, MD5 and date.
  Re-running the command skips what is complete.

A reference with ALT contigs and no ALT contig list stops the run: reads in the
regions those contigs duplicate would get mapping quality 0 and be lost.

## Run

From a project folder that contains `fastq/`:

```bash
bash hemaclone.sh run \
    --cores 32 --mem-gb 96 --batch-size 1 \
    --resources-dir /path/to/hg38_resources \
    --tmpdir /path/to/tmp \
    --set ch_filter.min_dp=100 --set ch_filter.min_vaf=0.02
```

- `--cores` and `--mem-gb` are hard limits on the container.
- `--batch-size N` runs N samples at a time; each batch finishes before the next starts.
- `--dry-run` shows the planned jobs; `bash hemaclone.sh run --help` lists every flag.
- Any setting in `config/config.yaml` can be changed with `--set KEY=VALUE`.
- The settings of every run are recorded in `logs/batches/`.

## Main outputs

| File | Content |
|---|---|
| `chip/chip_samples.tsv` | one row per sample: CHIP status, number of calls, largest VAF, genes |
| `chip/chip_calls.maf` | qualifying calls (maftools) |
| `chip/chip_calls.oncoplot.maf` | the same, plus a placeholder row per CHIP- sample for oncoplots |
| `chip/chip_candidates.tsv` | every candidate call with one column per criterion |
| `chip/chip_filter_counts.tsv` | per criterion: calls failing it, and calls failing only it |
| `qc/multiqc/chip_qc.html`, `qc/qc_flags.tsv` | QC report and per-sample flags |

## Report

`analysis/ch_analysis.Rmd` builds a maftools report from the outputs:

```bash
Rscript -e 'rmarkdown::render("analysis/ch_analysis.Rmd", params = list(project_dir = "/path/to/project"))'
```

## Interpretation

Without a matched normal, a call in blood can also come from tumour DNA,
therapy-selected clones, mosaicism, rare germline variants or sequencing
artefacts. Report calls as somatic blood variants consistent with CH, and
review recurrent positions, contamination and depth before biological
interpretation.

## Citation

Manuscript in preparation.

## License

MIT, see [`LICENSE`](LICENSE).
