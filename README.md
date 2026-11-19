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
- hg38 resources in one folder (GATK resource bundle names by default):
  reference FASTA with `.fai`, `.dict` and bwa-mem2 index; exome target
  interval list; CH gene panel BED; AF-only gnomAD; small ExAC common
  variants; dbSNP; Mills and 1000G indels; Funcotator data sources; optionally
  the 1000 Genomes panel of normals. Every file name can be changed in
  `config/config.yaml` or at launch.

## Build the image

```bash
git clone https://github.com/KevinDelGaudio/HemaClone.git
cd HemaClone
bash workflow/hemaclone.sh build          # builds hemaclone:<VERSION>
```

## Run

From a project folder that contains `fastq/`:

```bash
docker run --rm hemaclone:1.5.0 cat /opt/hemaclone/workflow/hemaclone.sh > hemaclone.sh

bash hemaclone.sh run \
    --cores 32 --mem-gb 96 --batch-size 1 \
    --resources-dir /path/to/hg38_resources \
    --tmpdir /path/to/tmp \
    --set reference.gene_panel_bed=chip_genes.exome_pad50.bed \
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
