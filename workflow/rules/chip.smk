# =============================================================================
# CHIP classification of the ch_final calls (criteria: chip block of
# config.yaml). Per sample, then one cohort set of tables and MAFs that
# maftools reads directly:
#   maf <- read.maf("chip/chip_calls.oncoplot.maf", clinicalData = "chip/chip_samples.tsv")
# =============================================================================

rule chip_sample:
    input:
        maf = "maf/{sample}_ch_final.funcotated.maf",
        vcf = "vcf/gnomad/{sample}.ch_final.vcf.gz",
    output:
        candidates = "chip/samples/{sample}.candidates.tsv",
        status = "chip/samples/{sample}.status.tsv",
    params:
        chip = config["chip"],
        min_dp = config["ch_filter"]["min_dp"],
        min_vaf = config["ch_filter"]["min_vaf"],
        pon_used = PON_ANNOTATE,
    threads: threads_of("chip_sample")
    resources:
        mem_mb = mem_of("chip_sample"),
    log: "logs/chip_sample/{sample}.log"
    script:
        "../scripts/chip_sample.py"


rule chip_cohort:
    input:
        tables = expand("chip/samples/{sample}.candidates.tsv", sample=SAMPLES),
        status = expand("chip/samples/{sample}.status.tsv", sample=SAMPLES),
    output:
        candidates = "chip/chip_candidates.tsv",
        maf = "chip/chip_calls.maf",
        oncoplot_maf = "chip/chip_calls.oncoplot.maf",
        samples = "chip/chip_samples.tsv",
    params:
        samples = SAMPLES,
        large_clone_vaf = config["chip"]["large_clone_vaf"],
    threads: threads_of("chip_cohort")
    resources:
        mem_mb = mem_of("chip_cohort"),
    log: "logs/chip_cohort/chip_cohort.log"
    script:
        "../scripts/chip_cohort.py"
