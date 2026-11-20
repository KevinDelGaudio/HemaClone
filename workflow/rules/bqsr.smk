# =============================================================================
# Base quality score recalibration, fit over the exome capture intervals (WES)
# or the WGS calling regions (WGS), so the error model has enough sites to be
# stable. In WGS the recalibrated BAM holds only the reads the later steps
# use: the gene panel (Mutect2) and the common sites (contamination), padded.
# =============================================================================


rule base_recalibrator:
    input:
        bam = "bam/{sample}.markdup.bam",
        bai = "bam/{sample}.markdup.bai",
    output:
        "bqsr/{sample}.recal.table",
    params:
        java = java_opts,
        ref = REF["fasta"],
        dbsnp = REF["dbsnp"],
        mills = REF["mills"],
        intervals = BQSR_INTERVALS,
        extra = TP["base_recalibrator"]["extra"],
    threads: threads_of("base_recalibrator")
    resources:
        mem_mb = mem_of("base_recalibrator"),
    log: "logs/base_recalibrator/{sample}.log"
    benchmark: "benchmarks/base_recalibrator/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' BaseRecalibrator"
        " -R {params.ref}"
        " -I {input.bam}"
        " --known-sites {params.dbsnp}"
        " --known-sites {params.mills}"
        " -L {params.intervals}"
        " {params.extra}"
        " -O {output}"
        " > {log} 2>&1"


rule apply_bqsr:
    input:
        bam = "bam/{sample}.markdup.bam",
        bai = "bam/{sample}.markdup.bai",
        table = "bqsr/{sample}.recal.table",
    output:
        bam = "bqsr/{sample}.bqsr.bam",
        bai = "bqsr/{sample}.bqsr.bai",
    params:
        java = java_opts,
        ref = REF["fasta"],
        intervals = apply_bqsr_intervals(),
        extra = TP["apply_bqsr"]["extra"],
    threads: threads_of("apply_bqsr")
    resources:
        mem_mb = mem_of("apply_bqsr"),
    log: "logs/apply_bqsr/{sample}.log"
    benchmark: "benchmarks/apply_bqsr/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' ApplyBQSR"
        " -R {params.ref}"
        " -I {input.bam}"
        " --bqsr-recal-file {input.table}"
        "{params.intervals}"
        " {params.extra}"
        " -O {output.bam}"
        " > {log} 2>&1"
