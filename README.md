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

### BAM or CRAM input

Put `<sample>.bam` or `<sample>.cram` files in the FASTQ folder (alone or next
to FASTQs). Their reads are extracted, with the original base qualities from
the `OQ` tag when present, and re-aligned with the same steps as FASTQ input,
so every sample is processed identically whatever reference the file was made
with. The sample name comes from the file name; the library from the header.
CRAM needs the FASTA it was written against:

```bash
bash hemaclone.sh run \
    --cores 32 --mem-gb 96 --batch-size 1 \
    --resources-dir /path/to/hg38_resources \
    --tmpdir /path/to/tmp \
    --set input.cram_reference=/path/to/cram_reference.fasta
```

Only read pairs are re-aligned: unpaired reads are dropped, and reads whose
mate is missing are counted in `logs/bam_remap/` as singletons. FastQC runs on
FASTQ input only.

### Sex and sample identity

somalier compares every pair of samples at known polymorphic sites:
`POSSIBLE_DUPLICATE` in `qc/qc_flags.tsv` marks pairs with relatedness of at
least 0.8 (the same individual, or a sample swap); `related_samples` lists
pairs from 0.4 (relatives). The X and Y depth ratios are reported for every
sample. Once male and female samples separate on a first run, set the Y ratio
cut-off and, optionally, a table of expected sex (sample, male/female, kept
with your sample lists) to flag `SEX_MISMATCH`:

```bash
    --set qc.sex_y_ratio_male=0.2 --set qc.expected_sex=/path/to/expected_sex.tsv
```

### Whole-genome data

Add `--data-type wgs` and choose the CH depth and VAF floors for your depth;
the run stops without them, because the defaults (DP 100, VAF 0.02) are set
for exome depth:

```bash
bash hemaclone.sh run \
    --cores 32 --mem-gb 96 --batch-size 1 \
    --resources-dir /path/to/hg38_resources \
    --tmpdir /path/to/tmp \
    --data-type wgs \
    --set ch_filter.min_dp=20 --set ch_filter.min_vaf=0.10
```

Probability of at least 3 alternate reads (the minimum for a call), by depth
and VAF:

| Depth | VAF 0.05 | VAF 0.10 | VAF 0.20 |
|---|---|---|---|
| 30x | 0.19 | 0.59 | 0.96 |
| 60x | 0.58 | 0.95 | 1.00 |
| 100x | 0.88 | 1.00 | 1.00 |

In WGS mode BQSR is fitted on the WGS calling regions, the recalibrated BAM
keeps only the gene panel and the contamination sites, the duplicate-marked
BAM is deleted once QC and BQSR are done (about 100 GB per 30x sample), and QC
uses genome-wide coverage (`qc.thresholds_wgs`). The exome interval list is
not needed.
- The settings of every run are recorded in `logs/batches/`.

## Main outputs

| File | Content |
|---|---|
| `chip/chip_samples.tsv` | one row per sample: CHIP status, number of calls, largest VAF, genes |
| `chip/chip_calls.maf` | qualifying calls (maftools) |
| `chip/chip_calls.oncoplot.maf` | the same, plus a placeholder row per CHIP- sample for oncoplots |
| `chip/chip_candidates.tsv` | every candidate call with one column per criterion |
| `chip/chip_filter_counts.tsv` | per criterion: calls failing it, and calls failing only it |
| `qc/multiqc/chip_qc.html`, `qc/qc_flags.tsv` | QC report and per-sample flags, sex and identity check included |
| `qc/somalier/relate.html` | somalier relatedness and sex report |
| `reports/ch_report.html` | CH report (maftools): methods with the thresholds applied, filtering, strand bias, indels in repeats, calls per sample, VAF, recurrent positions, oncoplot, lollipop plots |
| `reports/qc_report.html` | QC report: coverage against the CH depth floor, breadth, duplicates, insert size, contamination, sex and relatedness, input type |
| `reports/methods_report.html` | methods record: workflow version and command, tool versions, reference files with sources and MD5, every setting |
| `reports/results/*.csv` | tables of the reports |

## Report

Three reports are rendered at the end of a run, once every sample has
finished, with R from its own environment in the image: the CH report, the QC
report and the methods record. To render them again, for example after
excluding samples, run the same command with `--cohort-only`:

```bash
bash hemaclone.sh run \
    --cores 4 --mem-gb 16 --batch-size 1 \
    --resources-dir /path/to/hg38_resources \
    --tmpdir /path/to/tmp \
    --set ch_filter.min_dp=100 --set ch_filter.min_vaf=0.02 \
    --set report.exclude_samples=exclude.txt \
    --cohort-only
```

- Paths are read inside the container: relative paths start from the project
  folder; a file elsewhere needs its folder added with `--mount`.
- `report.exclude_samples`: a file of sample IDs, one per line (`#` starts a
  comment), left out of the CH report and its CSV tables and marked in the QC
  report. The tables in `chip/` keep every sample; the CH report names the
  samples it leaves out.
- `report.ch_rmd`, `report.qc_rmd`, `report.methods_rmd`: render your own
  version of a report. Copy the built-in one out of the image, edit it, and
  pass its path:

  ```bash
  docker run --rm hemaclone:1.6.0 cat /opt/hemaclone/workflow/report/ch_report.Rmd > my_ch_report.Rmd
  ```

  then add `--set report.ch_rmd=my_ch_report.Rmd` to the command.
  (`qc_report.Rmd` and `methods_report.Rmd` the same way).
- `report.enabled=false` skips the reports.
- If a report fails, the run says so and points to its log in `logs/`;
  every other output is complete.

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
