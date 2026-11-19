# =============================================================================
# QC: FastQC per FASTQ, samtools and Picard metrics per sample BAM,
# MultiQC report and per-sample threshold flags over the whole sheet.
# =============================================================================

rule fastqc:
    input:
        lambda wc: FASTQ_BY_BASE[wc.fq],
    output:
        html = "qc/fastqc/{fq}_fastqc.html",
        zip = "qc/fastqc/{fq}_fastqc.zip",
    params:
        run = QC_RUN,
    threads: threads_of("fastqc")
    resources:
        mem_mb = mem_of("fastqc"),
    log: "logs/fastqc/{fq}.log"
    benchmark: "benchmarks/fastqc/{fq}.tsv"
    shell:
        "{params.run} fastqc -t {threads} --outdir qc/fastqc --noextract --quiet {input}"
        " > {log} 2>&1"


rule samtools_qc:
    input:
        bam = "bam/{sample}.markdup.bam",
        bai = "bam/{sample}.markdup.bai",
    output:
        stats = "qc/bam/{sample}.stats.txt",
        flagstat = "qc/bam/{sample}.flagstat.txt",
        idxstats = "qc/bam/{sample}.idxstats.txt",
    threads: threads_of("samtools_qc")
    resources:
        mem_mb = mem_of("samtools_qc"),
    log: "logs/samtools_qc/{sample}.log"
    benchmark: "benchmarks/samtools_qc/{sample}.tsv"
    shell:
        "( samtools stats -@ {threads} {input.bam} > {output.stats}"
        " && samtools flagstat -@ {threads} {input.bam} > {output.flagstat}"
        " && samtools idxstats {input.bam} > {output.idxstats} ) 2> {log}"


rule alignment_summary:
    input:
        bam = "bam/{sample}.markdup.bam",
        bai = "bam/{sample}.markdup.bai",
    output:
        "qc/bam/{sample}.alignment_summary_metrics.txt",
    params:
        java = java_opts,
        ref = REF["fasta"],
    threads: threads_of("alignment_summary")
    resources:
        mem_mb = mem_of("alignment_summary"),
    log: "logs/alignment_summary/{sample}.log"
    benchmark: "benchmarks/alignment_summary/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' CollectAlignmentSummaryMetrics"
        " -I {input.bam}"
        " -R {params.ref}"
        " -O {output}"
        " --VALIDATION_STRINGENCY LENIENT"
        " > {log} 2>&1"


rule insert_size:
    input:
        bam = "bam/{sample}.markdup.bam",
        bai = "bam/{sample}.markdup.bai",
    output:
        metrics = "qc/bam/{sample}.insert_size_metrics.txt",
        histogram = "qc/bam/{sample}.insert_size_histogram.pdf",
    params:
        java = java_opts,
        min_pct = TP["insert_size"]["min_pct"],
    threads: threads_of("insert_size")
    resources:
        mem_mb = mem_of("insert_size"),
    log: "logs/insert_size/{sample}.log"
    benchmark: "benchmarks/insert_size/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' CollectInsertSizeMetrics"
        " -I {input.bam}"
        " -O {output.metrics}"
        " -H {output.histogram}"
        " -M {params.min_pct}"
        " --VALIDATION_STRINGENCY LENIENT"
        " > {log} 2>&1"


rule hs_metrics:
    input:
        bam = "bam/{sample}.markdup.bam",
        bai = "bam/{sample}.markdup.bai",
    output:
        "qc/bam/{sample}.hs_metrics.txt",
    params:
        java = java_opts,
        ref = REF["fasta"],
        intervals = REF["exome_intervals"],
        coverage_cap = TP["hs_metrics"]["coverage_cap"],
        min_mapq = TP["hs_metrics"]["min_mapq"],
        min_baseq = TP["hs_metrics"]["min_baseq"],
    threads: threads_of("hs_metrics")
    resources:
        mem_mb = mem_of("hs_metrics"),
    log: "logs/hs_metrics/{sample}.log"
    benchmark: "benchmarks/hs_metrics/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' CollectHsMetrics"
        " -I {input.bam}"
        " -R {params.ref}"
        " -BAIT_INTERVALS {params.intervals}"
        " -TARGET_INTERVALS {params.intervals}"
        " -O {output}"
        " -COVERAGE_CAP {params.coverage_cap}"
        " -MINIMUM_MAPPING_QUALITY {params.min_mapq}"
        " -MINIMUM_BASE_QUALITY {params.min_baseq}"
        " --VALIDATION_STRINGENCY LENIENT"
        " > {log} 2>&1"


rule multiqc:
    input:
        fastqc = expand("qc/fastqc/{fq}_fastqc.zip", fq=FASTQ_BY_BASE),
        samtools = expand("qc/bam/{sample}.stats.txt", sample=SAMPLES),
        alignment = expand("qc/bam/{sample}.alignment_summary_metrics.txt", sample=SAMPLES),
        insert_size = expand("qc/bam/{sample}.insert_size_metrics.txt", sample=SAMPLES),
        hs = expand("qc/bam/{sample}.hs_metrics.txt", sample=SAMPLES),
        dup = expand("qc/markdup/{sample}.dup_metrics.txt", sample=SAMPLES),
        config = MULTIQC_CONFIG,
    output:
        "qc/multiqc/chip_qc.html",
    params:
        run = QC_RUN,
    threads: threads_of("multiqc")
    resources:
        mem_mb = mem_of("multiqc"),
    log: "logs/multiqc/multiqc.log"
    benchmark: "benchmarks/multiqc/multiqc.tsv"
    shell:
        "{params.run} multiqc qc/fastqc qc/bam qc/markdup"
        " --config {input.config}"
        " --outdir qc/multiqc"
        " --filename chip_qc"
        " --interactive"
        " --force"
        " > {log} 2>&1"


rule qc_flags:
    input:
        hs = expand("qc/bam/{sample}.hs_metrics.txt", sample=SAMPLES),
        insert_size = expand("qc/bam/{sample}.insert_size_metrics.txt", sample=SAMPLES),
        dup = expand("qc/markdup/{sample}.dup_metrics.txt", sample=SAMPLES),
    output:
        "qc/qc_flags.tsv",
    params:
        samples = SAMPLES,
        thresholds = config["qc"]["thresholds"],
    threads: threads_of("qc_flags")
    resources:
        mem_mb = mem_of("qc_flags"),
    log: "logs/qc_flags/qc_flags.log"
    benchmark: "benchmarks/qc_flags/qc_flags.tsv"
    script:
        "../scripts/qc_flags.py"
