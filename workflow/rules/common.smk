# =============================================================================
# Sample sheet loading, preflight checks and helper functions.
# Everything here runs when Snakemake parses the workflow, including dry runs,
# so a broken sheet, a missing resource or an unset budget stops the run
# before any job starts.
# =============================================================================

import os
import re
import shutil
import subprocess
import sys

import pandas as pd


class PreflightError(Exception):
    pass


def log_preflight(text):
    print(f"[preflight] {text}", file=sys.stderr)


def fail_if(items, message):
    if items:
        listing = "\n  ".join(str(i) for i in items)
        raise PreflightError(f"{message}:\n  {listing}")


# -----------------------------------------------------------------------------
# Machine budget: must be set explicitly (profile or command line)
# -----------------------------------------------------------------------------
budget_mem = workflow.global_resources.get("mem_mb")
budget_mem = getattr(budget_mem, "value", budget_mem)
budget_cores = workflow.global_resources.get("_cores")
budget_cores = getattr(budget_cores, "value", budget_cores)

if budget_mem is None:
    raise PreflightError("No memory budget set. Launch with: bash workflow/hemaclone.sh run ... "
                         "(or, without Docker, python workflow/scripts/run_batches.py ...).")

oversized = [f"{rule}: {size['mem_mb']} MB" for rule, size in config["job_sizes"].items()
             if size["mem_mb"] > budget_mem]
fail_if(oversized, f"Jobs larger than the memory budget of {budget_mem} MB")

log_preflight(f"budget: {budget_cores} cores, {budget_mem} MB")


# -----------------------------------------------------------------------------
# Sample sheet
# -----------------------------------------------------------------------------
REQUIRED_COLUMNS = ["sample", "unit", "fq1", "fq2", "flowcell", "lane",
                    "barcode", "rg_id", "rg_pu", "lb"]
NAME_PATTERN = r"[A-Za-z0-9_-]+"
RG_PATTERN = r"[^\s'\"]+"

if not os.path.isfile(config["samples"]):
    raise PreflightError(f"Sample sheet not found: {config['samples']}. "
                         "Create it with workflow/scripts/make_samplesheet.py")

units = pd.read_csv(config["samples"], sep="\t", dtype=str, keep_default_na=False)

fail_if([c for c in REQUIRED_COLUMNS if c not in units.columns],
        "Sample sheet is missing columns")

empty_rows = units.index[units[REQUIRED_COLUMNS].eq("").any(axis=1)].tolist()
fail_if(empty_rows, "Sample sheet rows with empty cells (0-based)")

fail_if(sorted(units.loc[~units["sample"].str.fullmatch(NAME_PATTERN), "sample"].unique()),
        "Sample names with characters outside [A-Za-z0-9_-]")
fail_if(sorted(units.loc[~units["unit"].str.fullmatch(NAME_PATTERN), "unit"].unique()),
        "Unit names with characters outside [A-Za-z0-9_-]")

for column in ["rg_id", "rg_pu", "lb"]:
    fail_if(units.loc[~units[column].str.fullmatch(RG_PATTERN), column].tolist(),
            f"Read group field '{column}' contains whitespace or quotes")

duplicated_units = units.loc[units.duplicated(["sample", "unit"], keep=False), ["sample", "unit"]]
fail_if(duplicated_units.to_records(index=False).tolist(), "Duplicated (sample, unit) pairs")

duplicated_rg = units.loc[units.duplicated(["sample", "rg_id"], keep=False), ["sample", "rg_id"]]
fail_if(duplicated_rg.to_records(index=False).tolist(), "Read group IDs repeated within a sample")

all_fastq = pd.concat([units["fq1"], units["fq2"]], ignore_index=True)
fail_if(sorted(all_fastq[all_fastq.duplicated()].unique()), "FASTQ files listed more than once")
fail_if([p for p in all_fastq if not os.path.isfile(p)], "FASTQ files not found")

sample_names = pd.Series(sorted(units["sample"].unique()))
fail_if(sample_names[sample_names.str.lower().duplicated(keep=False)].tolist(),
        "Sample names that differ only by case")

libraries_per_sample = units.groupby("sample")["lb"].nunique()
for sample in libraries_per_sample[libraries_per_sample > 1].index:
    log_preflight(f"WARNING: sample {sample} has more than one library (LB); "
                  "duplicates are marked per library")

units = units.set_index(["sample", "unit"], drop=False).sort_index()
SAMPLES = sorted(units["sample"].unique())
UNIT_NAMES = sorted(units["unit"].unique())

# FastQC names its output after the FASTQ file name without extension
FASTQ_BY_BASE = {}
for path in all_fastq:
    base = re.sub(r"\.(fastq|fq)(\.gz)?$", "", os.path.basename(path))
    if base in FASTQ_BY_BASE:
        raise PreflightError(f"Two FASTQ files share the base name {base}: "
                             f"{FASTQ_BY_BASE[base]} and {path}")
    FASTQ_BY_BASE[base] = path

units_per_sample = units.groupby(level="sample").size()
log_preflight(f"sample sheet {config['samples']}: {len(units)} units, {len(SAMPLES)} samples, "
              f"{int((units_per_sample > 1).sum())} samples with more than one unit")


# -----------------------------------------------------------------------------
# Reference resources
# -----------------------------------------------------------------------------
REF_DIR = config["reference"]["dir"]
if not os.path.isdir(REF_DIR):
    raise PreflightError(f"Resources folder not found: {REF_DIR} (set it with --resources-dir)")

# File names are looked up inside reference.dir; absolute paths are used as is
REF = {}
for key, value in config["reference"].items():
    if key == "dir":
        continue
    REF[key] = value if os.path.isabs(value) else os.path.join(REF_DIR, value)

log_preflight(f"resources folder: {REF_DIR}")


def vcf_index(path):
    if path.endswith(".gz"):
        return f"{path}.tbi"
    return f"{path}.idx"


# Sequencing design: WES fits BQSR on the exome intervals and reports capture
# metrics; WGS uses the WGS calling regions and genome-wide coverage metrics.
DATA_TYPE = config["data_type"]
if DATA_TYPE not in ("wes", "wgs"):
    raise PreflightError(f"data_type must be 'wes' or 'wgs', got '{DATA_TYPE}'")
WGS = DATA_TYPE == "wgs"
BQSR_INTERVALS = REF["wgs_intervals"] if WGS else REF["exome_intervals"]
COVERAGE_METRICS = "wgs_metrics" if WGS else "hs_metrics"
log_preflight(f"data type: {DATA_TYPE}")

required_files = {
    "reference FASTA": REF["fasta"],
    "FASTA index": f"{REF['fasta']}.fai",
    "sequence dictionary": f"{os.path.splitext(REF['fasta'])[0]}.dict",
    "bwa-mem2 index": f"{REF['fasta']}.bwt.2bit.64",
    ("WGS intervals" if WGS else "exome intervals"): BQSR_INTERVALS,
    "gene panel BED": REF["gene_panel_bed"],
    "repeats BED": REF["repeats_bed"],
    "repeats BED index": f"{REF['repeats_bed']}.tbi",
}
for key in ["gnomad", "common_variants", "dbsnp", "mills"]:
    required_files[key] = REF[key]
    required_files[f"{key} index"] = vcf_index(REF[key])

missing_resources = [f"{label}: {path}" for label, path in required_files.items()
                     if not os.path.isfile(path)]
if not os.path.isdir(REF["funcotator_ds"]):
    missing_resources.append(f"Funcotator data sources: {REF['funcotator_ds']}")
fail_if(missing_resources, "Reference resources not found")

# Indel thresholds act on calls that already passed ch_filter.min_vaf and
# chip.min_alt_reads, so the higher of the two applies. A specific threshold
# below the general one (e.g. the 0.02 indel default under a WGS VAF floor of
# 0.10) has no effect: noted in the log, not an error.
CHIP = config["chip"]
for key, floor_key, floor in [
        ("indel_min_vaf", "ch_filter.min_vaf", config["ch_filter"]["min_vaf"]),
        ("repeat_indel_min_vaf", "ch_filter.min_vaf", config["ch_filter"]["min_vaf"]),
        ("indel_min_alt_reads", "chip.min_alt_reads", CHIP["min_alt_reads"]),
        ("repeat_indel_min_alt_reads", "chip.min_alt_reads", CHIP["min_alt_reads"])]:
    if CHIP[key] < floor:
        log_preflight(f"NOTE: chip.{key} = {CHIP[key]} is below {floor_key} = {floor}: "
                      f"{floor_key} applies to these indels")
if CHIP["max_sor"] <= 0:
    raise PreflightError(f"chip.max_sor must be positive, got {CHIP['max_sor']}")

# ALT and fix-patch contigs need bwa-mem2's ALT contig list (<FASTA>.alt): without
# it they count as ordinary chromosomes, reads in the regions they duplicate get
# mapping quality 0, and Mutect2 drops them.
alt_contigs = 0
with open(required_files["sequence dictionary"]) as handle:
    for line in handle:
        if line.startswith("@SQ"):
            name = line.split("SN:", 1)[1].split("\t", 1)[0].strip()
            if name.endswith(("_alt", "_fix")):
                alt_contigs += 1
alt_file = f"{REF['fasta']}.alt"
if alt_contigs and not os.path.isfile(alt_file):
    if not config["reference_checks"]["allow_alt_without_index"]:
        raise PreflightError(
            f"The reference has {alt_contigs} ALT or fix-patch contigs but no ALT contig list "
            f"({alt_file}). Without it, reads in the regions these contigs duplicate get "
            "mapping quality 0 and are lost before calling. Use the GATK reference built by "
            "'bash hemaclone.sh resources', or, to reproduce results made without the list, "
            "--set reference_checks.allow_alt_without_index=true")
    log_preflight(f"WARNING: {alt_contigs} ALT or fix-patch contigs without {alt_file}: "
                  "mapping quality 0 in the regions they duplicate (allowed by setting)")
elif alt_contigs:
    log_preflight(f"reference: {alt_contigs} ALT or fix-patch contigs, ALT contig list present")


# -----------------------------------------------------------------------------
# Panel of normals (optional, external)
# -----------------------------------------------------------------------------
SCRIPTS = os.path.join(WORKFLOW_HOME, "workflow", "scripts")

MULTIQC_CONFIG = config["qc"]["multiqc_config"]
if not os.path.isabs(MULTIQC_CONFIG):
    MULTIQC_CONFIG = os.path.join(WORKFLOW_HOME, MULTIQC_CONFIG)
fail_if([MULTIQC_CONFIG] if not os.path.isfile(MULTIQC_CONFIG) else [], "MultiQC config not found")

PON_VCF = config["pon"].get("vcf")
PON_MODE = config["pon"].get("mode", "annotate")

# Same lookup as the reference files: bare names inside reference.dir
if PON_VCF and not os.path.isabs(PON_VCF):
    PON_VCF = os.path.join(REF_DIR, PON_VCF)

if PON_VCF:
    if PON_MODE not in ("annotate", "filter"):
        raise PreflightError(f"pon.mode must be 'annotate' or 'filter', got '{PON_MODE}'")
    fail_if([p for p in [PON_VCF, f"{PON_VCF}.tbi"] if not os.path.isfile(p)],
            "Panel of normals files not found")
    log_preflight(f"panel of normals: {PON_VCF} (mode: {PON_MODE})")
else:
    log_preflight("panel of normals: none")

PON_FILTER = bool(PON_VCF) and PON_MODE == "filter"
PON_ANNOTATE = bool(PON_VCF) and PON_MODE == "annotate"


# -----------------------------------------------------------------------------
# Tools
# -----------------------------------------------------------------------------
TOOLS = config["tools"]
QC_ENV = config["qc"]["conda_env"]

needed_tools = [TOOLS["bwa"], "samtools", "gatk", "bcftools"]
if QC_ENV:
    needed_tools.append("conda")
else:
    needed_tools += ["fastqc", "multiqc"]
fail_if([t for t in needed_tools if shutil.which(t) is None], "Tools not found on PATH")

if QC_ENV:
    env_listing = subprocess.run(["conda", "env", "list"], capture_output=True, text=True).stdout
    env_names = [line.split()[0] for line in env_listing.splitlines()
                 if line.strip() and not line.startswith("#")]
    if QC_ENV not in env_names:
        raise PreflightError(f"Conda env '{QC_ENV}' (qc.conda_env) not found")
    QC_RUN = f"conda run --no-capture-output -n {QC_ENV}"
else:
    QC_RUN = ""

if TOOLS.get("gatk_version"):
    version_output = subprocess.run(["gatk", "--version"], capture_output=True, text=True)
    version_text = version_output.stdout + version_output.stderr
    if TOOLS["gatk_version"] not in version_text:
        raise PreflightError(f"GATK version mismatch: expected {TOOLS['gatk_version']}, got:\n{version_text}")


# -----------------------------------------------------------------------------
# Scratch space
# -----------------------------------------------------------------------------
os.makedirs(config["tmpdir"], exist_ok=True)
if not os.access(config["tmpdir"], os.W_OK):
    raise PreflightError(f"tmpdir is not writable: {config['tmpdir']}")

log_preflight(f"CH filters: DP >= {config['ch_filter']['min_dp']}, VAF >= {config['ch_filter']['min_vaf']}")


# -----------------------------------------------------------------------------
# Helpers used by the rules
# -----------------------------------------------------------------------------
TP = config["tool_params"]


def units_of(sample):
    return units.loc[sample, "unit"].tolist()


def unit_row(sample, unit):
    return units.loc[(sample, unit)]


def read_group(sample, unit):
    row = unit_row(sample, unit)
    return (f"@RG\\tID:{row['rg_id']}\\tSM:{sample}\\tPL:ILLUMINA"
            f"\\tLB:{row['lb']}\\tPU:{row['rg_pu']}")


def threads_of(rule_name):
    return config["job_sizes"][rule_name]["threads"]


def mem_of(rule_name):
    return config["job_sizes"][rule_name]["mem_mb"]


def java_opts(wildcards, resources):
    # Evaluated per job, so --set-resources overrides also change the heap
    heap_gb = max(1, int(resources.mem_mb * config["java_heap_fraction"] / 1024))
    return f"-Xmx{heap_gb}g -XX:+UseParallelGC -Djava.io.tmpdir={config['tmpdir']}"


def apply_bqsr_intervals():
    """WGS: ApplyBQSR writes only the gene panel and the contamination sites,
    padded; WES: the whole BAM (no -L)."""
    if not WGS:
        return ""
    padding = TP["apply_bqsr"]["wgs_padding"]
    return (f" -L {REF['gene_panel_bed']} -L {REF['common_variants']}"
            f" --interval-padding {padding}")


def filtered_vcf(wildcards):
    """FilterMutectCalls output, or its PoN-annotated copy when the PoN is in annotate mode."""
    if PON_ANNOTATE:
        return f"vcf/filtered/{wildcards.sample}.filtered.pon.vcf.gz"
    return f"vcf/filtered/{wildcards.sample}.filtered.vcf.gz"


def sample_outputs():
    outputs = []
    outputs += expand("maf/{sample}_all.funcotated.maf", sample=SAMPLES)
    outputs += expand("maf/{sample}_ch_final.funcotated.maf", sample=SAMPLES)
    outputs += expand("qc/fastqc/{fq}_fastqc.zip", fq=FASTQ_BY_BASE)
    for name in ["stats", "alignment_summary_metrics", "insert_size_metrics", COVERAGE_METRICS]:
        outputs += expand("qc/bam/{sample}.{name}.txt", sample=SAMPLES, name=name)
    outputs += expand("qc/markdup/{sample}.dup_metrics.txt", sample=SAMPLES)
    outputs += expand("chip/samples/{sample}.candidates.tsv", sample=SAMPLES)
    outputs += expand("chip/samples/{sample}.status.tsv", sample=SAMPLES)
    return outputs


def cohort_outputs():
    return ["qc/multiqc/chip_qc.html", "qc/qc_flags.tsv",
            "chip/chip_candidates.tsv", "chip/chip_calls.maf",
            "chip/chip_calls.oncoplot.maf", "chip/chip_samples.tsv",
            "chip/chip_filter_counts.tsv"]
