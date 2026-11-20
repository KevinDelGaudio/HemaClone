# =============================================================================
# Alignment: one bwa-mem2 job per unit (FASTQ pair), then MarkDuplicates
# merges all units of a sample into one BAM.
# =============================================================================

# Integrity of each FASTQ pair: both mates decompress, hold complete records,
# the same number of reads and the same read names in the same order.
# Not an input of bwa_map (adding one would re-align finished samples): the
# batch runner requests these files first, for unfinished samples only, and
# leaves a sample out of its batch when its check fails.
rule fastq_check:
    input:
        r1 = lambda wc: unit_row(wc.sample, wc.unit)["fq1"],
        r2 = lambda wc: unit_row(wc.sample, wc.unit)["fq2"],
    output:
        "qc/fastq_check/{sample}.{unit}.ok",
    params:
        script = os.path.join(SCRIPTS, "check_fastq_pair.sh"),
    threads: threads_of("fastq_check")
    resources:
        mem_mb = mem_of("fastq_check"),
    log: "logs/fastq_check/{sample}.{unit}.log"
    benchmark: "benchmarks/fastq_check/{sample}.{unit}.tsv"
    shell:
        "bash {params.script} {input.r1} {input.r2} {output} > {log} 2>&1"


rule bwa_map:
    input:
        r1 = lambda wc: unit_row(wc.sample, wc.unit)["fq1"],
        r2 = lambda wc: unit_row(wc.sample, wc.unit)["fq2"],
    output:
        temp("bam/units/{sample}.{unit}.sorted.bam"),
    params:
        bwa = TOOLS["bwa"],
        ref = REF["fasta"],
        rg = lambda wc: read_group(wc.sample, wc.unit),
        sort_mem = TP["bwa_map"]["sort_mem_per_thread"],
        sort_tmp = lambda wc: os.path.join(config["tmpdir"], f"{wc.sample}.{wc.unit}.sort"),
        extra = TP["bwa_map"]["extra"],
    threads: threads_of("bwa_map")
    resources:
        mem_mb = mem_of("bwa_map"),
    log: "logs/bwa_map/{sample}.{unit}.log"
    benchmark: "benchmarks/bwa_map/{sample}.{unit}.tsv"
    shell:
        "( {params.bwa} mem -t {threads} {params.extra} -R '{params.rg}' {params.ref} {input.r1} {input.r2}"
        " | samtools view -u -@ {threads} -"
        " | samtools sort -@ {threads} -m {params.sort_mem} -T {params.sort_tmp} -o {output} - )"
        " 2> {log}"


rule markdup:
    input:
        lambda wc: expand("bam/units/{sample}.{unit}.sorted.bam",
                          sample=wc.sample, unit=units_of(wc.sample)),
    # WGS: deleted once QC and BQSR have read it (about 100 GB per 30x sample)
    output:
        bam = temp("bam/{sample}.markdup.bam") if WGS else "bam/{sample}.markdup.bam",
        bai = temp("bam/{sample}.markdup.bai") if WGS else "bam/{sample}.markdup.bai",
        metrics = "qc/markdup/{sample}.dup_metrics.txt",
    params:
        inputs = lambda wc, input: " ".join(f"-I {bam}" for bam in input),
        java = java_opts,
        tmp = config["tmpdir"],
        extra = TP["markdup"]["extra"],
    threads: threads_of("markdup")
    resources:
        mem_mb = mem_of("markdup"),
    log: "logs/markdup/{sample}.log"
    benchmark: "benchmarks/markdup/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' MarkDuplicates"
        " {params.inputs}"
        " -O {output.bam}"
        " -M {output.metrics}"
        " --CREATE_INDEX true"
        " --TMP_DIR {params.tmp}"
        " {params.extra}"
        " > {log} 2>&1"
