# HemaClone: design and context

What HemaClone does, why it is built this way, the bugs already fixed, and
what is still open. Keep this file in sync with the code: update it in the
same commit as any change to steps, criteria, defaults or outputs.
Version history is in `CHANGELOG.md`.

---

## 1. Purpose

Reusable, Dockerized Snakemake workflow that calls clonal haematopoiesis (CH)
variants from blood exome sequencing without a matched normal, and classifies
each sample as CHIP+ or CHIP-. It replaces ad hoc bash scripts for mapping,
panel of normals, CH calling and QC with one reusable, versioned workflow.

Requirements set at design time:

- usable by the whole lab, through Docker; publishable, no lab or sensitive references
- every setting controllable from one complete launch command (no local settings file)
- strict control of how many samples run at once, and hard CPU and RAM limits
- FASTQ integrity checks; keep going past failed samples and list them
- lane-split and non-split FASTQs in the same run, without merging FASTQs
- auto-built sample sheet
- CHIP classification inside the workflow, outputs directly usable by maftools
- each sample classified from its own data only: results do not depend on which other samples are in the run

## 2. Repository layout

```
VERSION                       current version
CHANGELOG.md
CONTRIBUTING.md               rules for changing code, config and docs
docs/CONTEXT.md               this file
docs/COMPARATORS.md           existing CH tools for the comparison analysis
analysis/ch_analysis.Rmd      downstream maftools report
analysis/exclude_samples.sh   removes samples from the cohort tables
config/config.yaml            defaults for every setting
config/profile/config.yaml    Snakemake behaviour (keep-going, rerun-incomplete, latency-wait 60)
config/multiqc_config.yaml
workflow/Snakefile
workflow/hemaclone.sh         launcher (Docker): run | samplesheet | shell | build
workflow/docker/Dockerfile    micromamba base, code baked into /opt/hemaclone
workflow/docker/environment.yaml
workflow/rules/               common, mapping, qc, bqsr, calling, annotation, chip (.smk)
workflow/scripts/             make_samplesheet.py, run_batches.py, check_fastq_pair.sh,
                              qc_flags.py, chip_sample.py, chip_cohort.py
tests/make_fixture.sh         synthetic project (tiny FASTQs, empty resources) for dry runs
tests/resources_test.sh       builds a resources folder from a tiny bundle inside the image (CI)
tests/ci_checks.sh            checks run by CI against an image
tests/compare_outputs.sh      compares the results of two runs (outputs-identical changes)
.github/workflows/ci.yml      CI on pull requests and pushes to main and release branches
.github/workflows/release.yml image published to GHCR when a release is published
```

Project folder layout (where the workflow is launched): `fastq/` with the
input FASTQs; outputs split by type (`bam/`, `bqsr/`, `vcf/`, `orientation/`, `maf/`, `qc/`,
`contamination/`, `chip/`, `logs/`, `benchmarks/`).

## 3. Infrastructure

- Docker image `hemaclone:<VERSION>`; the tools layer is cached, so a code change does not reinstall the tools
- Pinned tools: Snakemake 9.27.1, GATK 4.6.2.0, bwa-mem2 2.2.1, samtools and bcftools 1.24, FastQC 0.12.1 (0.13 dropped `--noextract`), MultiQC 1.35, pandas 3.1.0, Python 3.12, R 4.4.1
- The run stops if `gatk --version` differs from `tools.gatk_version`
- `--dev` runs the code of a local checkout with the tools of the same image
- Get the launcher for a new project from the image:
  `docker run --rm hemaclone:<version> cat /opt/hemaclone/workflow/hemaclone.sh > hemaclone.sh`
- Build: `bash workflow/hemaclone.sh build` (adds the proxy build args when needed)
- Published images: each GitHub release `vX.Y.Z` builds the image and pushes it to `ghcr.io/kevindelgaudio/hemaclone:X.Y.Z`; the release stops if the tag is not `v` + `VERSION`. The launcher expects the local name, so after pulling: `docker tag ghcr.io/kevindelgaudio/hemaclone:X.Y.Z hemaclone:X.Y.Z`
- CI: every pull request builds the image (tools layer cached) and runs `tests/ci_checks.sh`: script syntax, launcher dry run and whole-graph dry run on a synthetic project, Snakemake lint (blocking; the per-rule container lint is ignored, since the whole workflow runs in one image). Nothing is aligned or called, so CI catches broken rules, settings and launcher logic, not wrong results

## 4. Launching

Settings precedence: `config.yaml` defaults < `--settings FILE` < `--set KEY=VALUE` < dedicated flags.
Every run writes its resolved settings to `logs/batches/run_<time>_settings.yaml`.

Main flags: `--cores`, `--mem-gb` (hard container limits), `--batch-size`,
`--only`, `--sample-sheet`, `--fastq-dir`, `--rebuild-sheet`,
`--allow-unparsed-headers`, `--resources-dir`, `--tmpdir`, `--settings`,
`--set`, `--dry-run`, `--cohort-only`, `--image`, `--dev`, `--mount`.

- Batches of N samples run to completion one after another (`--batch-size 1` = one sample at a time)
- Jobs are scheduled into `machine.schedule_fraction` (0.9) of `--mem-gb`, leaving headroom for Java overhead
- Java heap = `java_heap_fraction` (0.85) of each job's memory; `--set job_sizes.<rule>.mem_mb=N` resizes a job
- Before the first batch, the preflight checks run once on the full sample sheet; if they fail, the run stops with their message and no batch starts (a missing resource or bad setting concerns the run, not a sample)
- Failed samples do not stop the run; they are listed in `logs/batches/failed_<time>.tsv`
- Every rule writes a log (`logs/<rule>/`), including the error that stopped it, and a runtime and memory benchmark (`benchmarks/<rule>/`)

Example used in production (one sample per batch, 96 GB):

```bash
bash hemaclone.sh run --cores 32 --mem-gb 96 --batch-size 1 \
    --resources-dir /path/to/gatk_resources --tmpdir /path/to/tmp \
    --set reference.gene_panel_bed=chip_genes.exome_pad50.bed \
    --set job_sizes.bwa_map.mem_mb=32000
```

## 5. Steps

0. **Resources** (`bash hemaclone.sh resources`, `workflow/scripts/get_resources.py`), once per resources folder: GATK hg38 reference bundle (FASTA, `.fai`, `.dict`, ALT contig list saved as `<FASTA>.alt` for alt-aware bwa-mem2 mapping, dbSNP138, Mills/1000G indels, WGS calling regions) and somatic resources (gnomAD AF-only, small ExAC common sites, 1000 Genomes PoN), each checked against the server's size and MD5; bwa-mem2 index; Funcotator somatic data sources; GIAB tandem repeats and homopolymers BED (kept as downloaded under `giab/`, bgzipped and tabix-indexed as `reference.repeats_bed`); the kit BED (`--targets`) as an interval list; the CH gene panel BED = gene spans of `resources.panel_genes` from the Funcotator GENCODE GTF, intersected with the kit targets, each piece padded by `resources.panel_padding` (50 bp) and merged per gene. Genes missing from the GTF or found at two loci stop the step; genes without kit targets are listed. `MANIFEST.tsv` records source, size, MD5 and date of every file; finished files are skipped on re-runs.
1. **Sample sheet** built from `fastq/`. Lane-split FASTQs stay separate units with their own read groups; one BAM per sample. Preflight checks: required columns, name characters, duplicated units/read groups/FASTQs, case-only name clashes, missing resources, a reference with ALT or fix-patch contigs but no `<FASTA>.alt` (stops unless `reference_checks.allow_alt_without_index`), PoN files, tools, writable tmpdir, job sizes within the memory budget.
   BAM/CRAM files in the same folder are one unit each (`B1`, path in the sheet's `bam` column, empty `fq1`/`fq2`); sample name from the file name (a different header `SM` is a warning), `LB` from the header when it has one library. CRAM needs `input.cram_reference` (checked at preflight). Sheets without a `bam` column are FASTQ-only.
2. **FASTQ check** per pair: complete gzip decompression, matching read counts and read names in both mates. Truncated pairs fail the sample. BAM/CRAM units: `samtools quickcheck -u` (`bam_check`; unaligned BAMs accepted).
3. **Mapping**: bwa-mem2 per unit, samtools sort. BAM/CRAM units (`bam_remap`), streamed with no intermediate FASTQ: `samtools collate` | `samtools fastq -O` (original qualities from `OQ`; secondary, supplementary and unpaired reads dropped; reads whose mate is missing are counted in the log as singletons) | `bwa-mem2 mem -p` | sort; a file with no read pairs fails its sample. GATK MarkDuplicates per library after merging units. With `data_type: wgs` the duplicate-marked BAM is temporary: deleted once QC and BQSR have read it.
4. **BQSR**: BaseRecalibrator with dbSNP138 and Mills/1000G indels as known sites, over the exome targets (WES) or `reference.wgs_intervals` (WGS); ApplyBQSR. In WGS, ApplyBQSR writes only the reads over the gene panel and the `small_exac_common_3` sites, padded by `tool_params.apply_bqsr.wgs_padding` (150 bp): the only regions Mutect2 and the contamination estimate read.
5. **Calling**: Mutect2 tumor-only on the CHIP gene panel BED, gnomAD AF-only as germline resource, F1R2 collection.
6. **Artefact models**: LearnReadOrientationModel; GetPileupSummaries and CalculateContamination on `small_exac_common_3`; FilterMutectCalls with orientation priors and contamination table.
7. **Panel of normals**: `1000g_pon.hg38.vcf.gz` by default. `mode: annotate` adds an INFO/PON flag with bcftools and keeps the calls; `mode: filter` passes the PoN to Mutect2. Never use a PoN built from the cohort being analysed: recurrent CHIP hotspots would be treated as artefacts. Off for one run: `--set pon.vcf=null`.
8. **PASS only**: SelectVariants `--exclude-filtered`.
9. **CH genotype filters**: VariantFiltration `DP < min_dp` (LowDP) and `AF < min_vaf` (LowVAF), then SelectVariants `--exclude-filtered --set-filtered-gt-to-nocall --exclude-non-variants`.
10. **gnomAD annotation**: bcftools `annotate -c INFO/GNOMAD_AF:=INFO/AF` from the local AF-only gnomAD, allele-matched. Funcotator's own gnomAD sources are not used (they ship as unextracted tarballs).
11. **Funcotator** to MAF (hg38) for two tiers: `_all` (PASS) and `_ch_final` (after the CH filters).
12. **QC**: FastQC, samtools stats, CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, CollectHsMetrics (WES; coverage cap 1000, MAPQ and BQ at least 20) or CollectWgsMetrics over `reference.wgs_intervals` (WGS; coverage cap 250), somalier sex and identity check (`somalier extract` per sample from the duplicate-marked BAM, inside its batch; `somalier relate` over the cohort, `qc/somalier/relate.*`), MultiQC, `qc/qc_flags.tsv`.
13. **CHIP classification**: `chip_sample.py` per sample, `chip_cohort.py` for the cohort tables.

### Calling intervals

Two BED versions, chosen at launch, with no intervals logic hard-coded in the workflow:

- `chip_genes.exome_pad50.bed`: gene bodies of the CH panel intersected with the exome targets, padded by 50 bp (WES)
- `chip_genes.bed`: gene bodies (for a future WGS use)

## 6. CHIP classification criteria

All criteria are configurable at launch (`--set chip.<key>=...`, `--set ch_filter.<key>=...`).

| Criterion | Default |
|---|---|
| Protein-altering class | Missense, Nonsense, Nonstop, Frame_Shift_Del/Ins, In_Frame_Del/Ins, Splice_Site, Translation_Start_Site |
| Alternate reads | at least 3 |
| Depth | at least `ch_filter.min_dp` (100) |
| VAF floor | at least `ch_filter.min_vaf` (0.02) |
| VAF ceiling | at most `chip.max_vaf` (0.35) |
| Germline exclusion | one-sided binomial test VAF < 0.5, p < 0.01 |
| gnomAD AF | at most 0.01, variants absent from gnomAD = 0 |
| PoN | excluded if in the PoN, unless at a hotspot |
| Strand bias | strand odds ratio at most `chip.max_sor` (3), from Mutect2 `FORMAT/SB` (ref and alt reads on each strand), GATK formula with 1 added to each count |
| Indels | VAF at least `chip.indel_min_vaf` and alt reads at least `chip.indel_min_alt_reads`; defaults equal to the general ones (0.02, 3), so no effect unless raised |
| Indels in repeats | an indel flagged `STR` by Mutect2 or overlapping the GIAB tandem-repeat and homopolymer BED (`reference.repeats_bed`) needs VAF at least `chip.repeat_indel_min_vaf` (0.05) and alt reads at least `chip.repeat_indel_min_alt_reads` (5); no gene or hotspot is exempt (ASXL1 G646 included) |

PoN-exempt hotspots: DNMT3A R882, JAK2 V617F, SF3B1 K700/K666, SRSF2 P95,
U2AF1 S34/Q157, IDH1 R132, IDH2 R140/R172, GNB1 K57, KRAS and NRAS G12/G13/Q61.

- A sample is **CHIP+** with at least one qualifying call in a gene not in `status_exclude_genes` (IL1A, IL1B, IL18, IL33, NLRP3, P2RX7: reported, not counted).
- **Large clone**: VAF at least 0.10.
- MAF rows are matched one by one to VCF alleles (Funcotator writes one row per ALT allele in VCF order); a mismatch stops the step.
- Indel and repeat-indel thresholds act on calls that already passed `ch_filter.min_vaf` and `chip.min_alt_reads`, so the higher of the two applies; a specific threshold below the general one has no effect, noted in the preflight log (e.g. the indel defaults under a WGS VAF floor of 0.10).
- Strand counts are one set per VCF site, shared by the alleles of a multi-allelic site. A VCF without `FORMAT/SB` stops the step; without the `INFO/STR` header line, repeats come from the GIAB BED only (logged).
- An indel's span for the BED lookup: the anchor base, the deleted bases and the next base (where an insertion goes).
- These filters add to the Mutect2 strand and polymerase-slippage filters that PASS calls already passed.
- gnomAD availability is read from the VCF header. A sample without calls has no INFO columns in its MAF, which is expected; a MAF with calls but no gnomAD column stops the step.

## 7. Outputs

| File | Content |
|---|---|
| `chip/chip_samples.tsv` | one row per sample, CHIP- included: status, n calls, max and summed VAF, large clone, genes, gnomad_checked |
| `chip/chip_calls.maf` | counted calls only, maftools-ready |
| `chip/chip_calls.oncoplot.maf` | the same plus one Silent `CHIP_NEGATIVE_PLACEHOLDER` row per CHIP- sample, so oncoplot shows the whole cohort |
| `chip/chip_candidates.tsv` | every ch_final call with one `pass_` column per criterion, `qualifies`, `counts_for_status`, and the values behind the new criteria: `sb_ref_fwd`, `sb_ref_rev`, `sb_alt_fwd`, `sb_alt_rev`, `sor`, `is_indel`, `in_str`, `in_giab_repeat`, `in_repeat` |
| `chip/chip_filter_counts.tsv` | per criterion: candidate calls, calls failing it, calls failing it and no other (those it alone removes) |
| `chip/samples/<s>.candidates.tsv`, `<s>.status.tsv` | per-sample tables |
| `contamination/<s>.contamination.table` | CalculateContamination estimate |
| `qc/multiqc/chip_qc.html`, `qc/qc_flags.tsv` | QC report and per-sample flags |

WGS QC flags (`qc.thresholds_wgs`, columns `mean_coverage`, `pct_20x`, `pct_exc_dupe`, `pct_exc_mapq`, `median_insert_size`, `pct_dup_markdup`): mean coverage at least 25, at least 80% of the calling regions at 20x, duplication at most 30%.

Sex and identity columns (WES and WGS): `x_depth_ratio`, `y_depth_ratio` (mean depth at the somalier sites on X and Y over the autosomal mean), `inferred_sex` (male when `y_depth_ratio` is at least `qc.sex_y_ratio_male`; `NA` while that is null), `related_samples` (other samples with relatedness at least `qc.min_relatedness_list`, 0.4). Flags: `POSSIBLE_DUPLICATE` (relatedness at least `qc.min_relatedness_flag`, 0.8: the same individual or a sample swap) and `SEX_MISMATCH` (inferred sex differs from `qc.expected_sex`, an optional run-time table). somalier 0.3.1 runs from its own environment in the image (`workflow/docker/somalier.yaml`), because its htslib pin differs from that of samtools and bcftools.

WES QC flag thresholds: mean coverage at least 30, at least 80% of target at 20x,
at least 40% selected bases, duplication at most 50%, dropout at most 10.

## 8. Downstream analysis: `ch_analysis.Rmd`

maftools-only notebook, no clinical or transcriptomic data required.

- Settings and thresholds are read from the latest non-dry-run `logs/batches/run_*_settings.yaml` with a chip block, so the methods text reports the values actually applied.
- Stops if any sample has `gnomad_checked = FALSE`, or if the samples differ between the oncoplot MAF and the sample table.
- Sections: methods (steps and criteria, inline values), filtering overview (calls failing each criterion, and failing only that one), VAF distribution, recurrent positions, prevalence, plotmafSummary, oncoplot (`annotationOrder = c("CHIP+", "CHIP-")` puts negatives on the right; `barcode_mar = 10`), plotVaf, lollipop plots, CSV export to `results/`.
- maftools cannot load a MAF made of placeholders only, so MAFs are loaded only when at least one sample is CHIP+.
- Lollipop plots: labels only for residues mutated in at least two samples; leaving out `labelPos` removes all labels.
- `params$project_dir` must point to the folder that holds `chip/` and `logs/`; it becomes `knitr::opts_knit$set(root.dir = ...)`.
- Render with the project folder as parameter:
  `Rscript -e 'rmarkdown::render("ch_analysis.Rmd", params = list(project_dir = "/path/to/project"))'`

Companion script `exclude_samples.sh --exclude IDS_FILE PROJECT_DIR`: removes
the samples listed in IDS_FILE (one per line, kept outside the repository)
from the four cohort tables, keeping the originals in `chip_unfiltered/`;
rerunnable, stops on unknown IDs.

## 9. Bugs fixed (details in CHANGELOG.md)

- Snakemake reserved rule name `insert` renamed `insert_size`.
- FastQC pinned to 0.12.1 (0.13 dropped `--noextract`).
- **CH filter was a no-op** in the original bash pipeline: VariantFiltration only set genotype FT, so the ch_final VCFs were identical to the PASS VCFs. Fixed with `--set-filtered-gt-to-nocall --exclude-non-variants`.
- Failed-sample list renamed `failed_<time>.tsv` (clashed with the `run_*.tsv` glob).
- `LC_ALL=C` in the FASTQ check (locales with comma decimals broke printf).
- 1.4.0: gnomAD criterion was never applied (Funcotator gnomAD tarballs unextracted) and `gnomad_checked` was vacuous. Replaced by bcftools `GNOMAD_AF` from the local AF-only gnomAD, plus a per-sample status file.
- 1.4.1: samples with zero candidates reported `gnomad_checked = False`. Now read from the VCF header; stop if calls exist but the column is missing.
- Paired-read name mismatch detected by the FASTQ check (bwa error "paired reads have different names") rather than failing inside bwa.

## 10. Known limitations and open items

- **Sensitivity is unmeasured.** Candidate benchmark: GIAB mixtures or dilutions with known low-VAF variants, sensitivity and false positives by VAF bin.
- **min_dp vs real depth**: check that `ch_filter.min_dp` fits the achieved panel depth (CollectHsMetrics) before trusting calls near the floor.
- **No sample identity check** (e.g. Conpair) between samples of the same individual.
- **Systematic artefacts pass the filters.** Recurrent indels in low-complexity regions of long genes (KMT2A, KMT2D, BCORL1, CREBBP, SRCAP, ZNF318) and missense calls in GC-rich genes (CEBPA, GATA2, WT1) can clear every criterion, because an external PoN does not contain artefacts specific to one lab or protocol. Signs: many calls concentrated in a few samples, the same genes and positions across unrelated individuals, VAF near the floor, frameshift and in-frame insertions. Candidate filters, discussed but not implemented:
  - indels in tandem repeats and homopolymers, and indels in general: stricter thresholds, implemented in 1.6.0 (section 6)
  - within-cohort recurrence flag for non-hotspot positions: not adopted, because each sample is classified from its own data only (section 1); recurrent positions stay a review item in `ch_analysis.Rmd`
- **Optional matched-tumor annotation**, discussed but not implemented: genotype each blood call in the same patient's tumor BAM (Mutect2 `--alleles` or pileup) and compare VAFs. Tumor VAF well above blood VAF suggests tumor-derived; lower or absent suggests CH. Annotation only, never automatic removal; the tumor must not be used as a matched normal, because CH infiltrating the tumor would be removed.
- IGV review of a few calls per suspicious gene before reporting.

## 11. Interpretation notes

- Tumor-only calling cannot fully separate CH from germline, mosaicism, tumor-derived DNA or artefacts. Report calls as somatic blood variants consistent with CH.
- The 0.35 VAF ceiling and the binomial test exclude heterozygous germline variants by design, so inherited predisposition variants are not assessed.
- In young or non-elderly cohorts CH is rare; a high CHIP+ rate or dozens of calls per sample point to a technical cause first.
- Tumor cells in blood DNA at the WES detection limit (about 4% of cells for a heterozygous call at 2% VAF) require extensive marrow involvement.

## 12. Version 1.6.0: scope and decisions

Changes go through `release/1.6.0`, one pull request per change; `VERSION`
becomes 1.6.0 when the release branch is merged into `main`.

| Area | Decision |
|---|---|
| Release process | release branch; pull requests into it do not bump `VERSION` |
| Infrastructure | `pyyaml` pinned; CI on pull requests; image published to GHCR on release; launcher keeps the local name `hemaclone:X.Y.Z` |
| Polish | conservative clean-up, outputs identical, checked on a test run |
| Resources | `hemaclone.sh resources` downloads and indexes the hg38 resources and writes a manifest with sources and checksums; the CH gene panel BED is built from a gene list kept in the repository; exome targets stay user-supplied |
| Reference | GATK `Homo_sapiens_assembly38` with its ALT contig list (alt-aware bwa-mem2 mapping); a reference with ALT or fix-patch contigs and no list stops the run |
| Gene panel | 81 genes in `resources.panel_genes`; gene spans from the Funcotator GENCODE GTF intersected with the kit targets, padded 50 bp |
| Strand bias | new criterion: strand odds ratio (SOR, from `FORMAT/SB`) at most `chip.max_sor` (3); reads per strand reported as columns, not filtered |
| Indels in repeats | stricter VAF and alt-read thresholds for indels in tandem repeats or homopolymers (Mutect2 `STR` and the GIAB tandem-repeat and homopolymer BED), not exclusion; defaults VAF 0.05 and 5 alt reads, to be tuned with the benchmark; ASXL1 G646 is subject to them, not exempt |
| Indels | separate VAF and alt-read thresholds for indels and SNVs |
| Within-cohort recurrence | not implemented: samples are classified independently |
| Low mappability | not implemented |
| Gene-specific variant allow-list | deferred to a later version |
| WGS | `data_type: wes` or `wgs` (flag `--data-type`) switches BQSR intervals and the QC metrics; the recalibrated BAM covers only the panel and the contamination sites; the duplicate-marked BAM is deleted after QC and BQSR; a WGS run must set `ch_filter.min_dp` and `ch_filter.min_vaf` at launch |
| Input | BAM or CRAM accepted and re-aligned (reads extracted with original qualities, mapped like FASTQ); sample name from the file name; FastQC on FASTQ input only |
| Sample QC | sex and sample identity check with somalier 0.3.1, in its own environment in the image; duplicate flag at relatedness 0.8; sex inferred only once `qc.sex_y_ratio_male` is set from a first cohort |
| Benchmark | after 1.6.0, against the tools in `docs/COMPARATORS.md` |
