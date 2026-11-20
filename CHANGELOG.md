# Changelog

All notable changes to HemaClone. Versions follow `MAJOR.MINOR.PATCH`: a
change to calls or classification bumps at least MINOR; a change that leaves
every output identical is a PATCH. The version in `VERSION`, the git tag
`vX.Y.Z` and the Docker image tag `hemaclone:X.Y.Z` always match.

## [Unreleased]

Changes merged into `release/1.6.0`.

### Added
- CI on pull requests and pushes to `main` and release branches
  (`.github/workflows/ci.yml`): builds the image and runs `tests/ci_checks.sh`
  (script syntax, launcher dry run and whole-graph dry run on a synthetic
  project from `tests/make_fixture.sh`, Snakemake lint reported).
- Image published to `ghcr.io/kevindelgaudio/hemaclone:<VERSION>` when a GitHub
  release is published (`.github/workflows/release.yml`); the release stops if
  its tag is not `v<VERSION>`.
- Benchmarks for `pon_annotate`, `gnomad_annotate`, `multiqc`, `qc_flags`,
  `chip_sample` and `chip_cohort`: every rule now records runtime and memory.
- `tests/compare_outputs.sh OLD_PROJECT NEW_PROJECT`: compares the VCF, MAF,
  CHIP and QC-flag records of two runs, skipping header lines with dates and
  versions.

### Changed
- `pyyaml` pinned in `workflow/docker/environment.yaml`, to the version already
  in the 1.5.0 image; environment renamed `hemaclone`. Tools unchanged.
- `CONTRIBUTING.md`: release-branch process; within-cohort recurrence described
  as a design choice (each sample classified independently), not a gap.
- Code clean-up, outputs identical (same planned jobs and shell commands; same
  classification tables on a test case):
  - `chip_sample.py`, `chip_cohort.py`, `qc_flags.py` write the error that
    stops them to the rule log (Snakemake does not capture a script's stderr
    there), and close their logs on every exit.
  - Environment variable `CHIP_MODE` renamed `HEMACLONE_MODE` (the recorded
    `workflow.mode` setting is unchanged).
  - `filtered_vcf()` moved to `common.smk` with the other helpers; file names
    built with f-strings.
  - `markdup` builds its `-I` arguments from its inputs.
  - Launcher help printed up to the end of the header comment instead of fixed
    line numbers.
  - `run_batches.py` runs the preflight checks once on the full sample sheet
    before the first batch and stops the run if they fail. Before, a missing
    resource marked every sample as failed with a misleading "FASTQ check
    failed" reason, batch after batch.
- Snakemake lint is blocking in CI (the per-rule container lint is ignored).

## [1.5.0]

Renamed to HemaClone. No change to calling, filtering or classification:
outputs are identical to 1.4.1.

### Changed
- Launcher `workflow/chip.sh` renamed `workflow/hemaclone.sh`; image
  `chip-workflow:<version>` renamed `hemaclone:<version>`; code inside the
  image moved from `/opt/chip-workflow` to `/opt/hemaclone`; environment
  variable `CHIP_IMAGE` renamed `HEMACLONE_IMAGE`.
- Output folder `chip/`, the `chip:` config block and the `CHIP_status`
  columns keep their names (CHIP is the clinical term they describe).

### Added
- `README.md`, `CONTRIBUTING.md`, `docs/CONTEXT.md`, `docs/COMPARATORS.md`, this changelog.
- `analysis/ch_analysis.Rmd`: maftools report; the project folder is a render
  parameter (`params$project_dir`).
- `analysis/exclude_samples.sh`: removes samples listed in a file
  (`--exclude IDS_FILE`) from the cohort tables, keeping the originals.

## [1.4.1]

### Fixed
- Samples with zero candidate calls reported `gnomad_checked = False`, because
  Funcotator writes no INFO columns for an empty VCF. gnomAD availability is
  now read from the VCF header; the step stops if calls exist but the gnomAD
  column is missing.

## [1.4.0]

### Fixed
- The gnomAD criterion was never applied: Funcotator's gnomAD data sources ship
  as unextracted tarballs, and `gnomad_checked` was true even without a gnomAD
  column. gnomAD AF is now added with bcftools from the local AF-only gnomAD
  VCF (`INFO/GNOMAD_AF`, allele-matched), with a per-sample status file.

### Added
- CHIP classification inside the workflow (`chip_sample.py`, `chip_cohort.py`):
  per-criterion pass columns, cohort tables, maftools-ready MAFs with Silent
  placeholders for CHIP- samples.

## Before 1.4.0

### Fixed
- The CH depth and VAF filter was a no-op: VariantFiltration only set the
  genotype FT field, so `ch_final` VCFs were identical to the PASS VCFs. Fixed
  with `SelectVariants --set-filtered-gt-to-nocall --exclude-non-variants`.
- Rule `insert` renamed `insert_size` (reserved name in Snakemake).
- FastQC pinned to 0.12.1 (0.13 dropped `--noextract`).
- Failed-sample list renamed `failed_<time>.tsv` (clashed with the `run_*.tsv` glob).
- `LC_ALL=C` in the FASTQ check (locales with comma decimals broke printf).
- Paired-read name mismatches are caught by the FASTQ check instead of failing inside bwa-mem2.
