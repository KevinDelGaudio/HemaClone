# =============================================================================
# Tumor-only Mutect2 on the CHIP gene panel, orientation and contamination
# models, FilterMutectCalls, PASS selection and CH-specific DP/VAF filters.
# =============================================================================

rule mutect2:
    input:
        bam = "bqsr/{sample}.bqsr.bam",
        bai = "bqsr/{sample}.bqsr.bai",
        pon = [PON_VCF] if PON_FILTER else [],
    output:
        vcf = "vcf/raw/{sample}.raw.vcf.gz",
        tbi = "vcf/raw/{sample}.raw.vcf.gz.tbi",
        stats = "vcf/raw/{sample}.raw.vcf.gz.stats",
        f1r2 = "orientation/{sample}.f1r2.tar.gz",
    params:
        java = java_opts,
        ref = REF["fasta"],
        gnomad = REF["gnomad"],
        panel = REF["gene_panel_bed"],
        pon_arg = f"--panel-of-normals {PON_VCF}" if PON_FILTER else "",
        extra = TP["mutect2"]["extra"],
    threads: threads_of("mutect2")
    resources:
        mem_mb = mem_of("mutect2"),
    log: "logs/mutect2/{sample}.log"
    benchmark: "benchmarks/mutect2/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' Mutect2"
        " -R {params.ref}"
        " -I {input.bam}"
        " --tumor-sample {wildcards.sample}"
        " --germline-resource {params.gnomad}"
        " --f1r2-tar-gz {output.f1r2}"
        " -L {params.panel}"
        " --native-pair-hmm-threads {threads}"
        " {params.pon_arg}"
        " {params.extra}"
        " -O {output.vcf}"
        " > {log} 2>&1"


rule learn_orientation:
    input:
        "orientation/{sample}.f1r2.tar.gz",
    output:
        "orientation/{sample}.read_orientation_model.tar.gz",
    params:
        java = java_opts,
    threads: threads_of("learn_orientation")
    resources:
        mem_mb = mem_of("learn_orientation"),
    log: "logs/learn_orientation/{sample}.log"
    benchmark: "benchmarks/learn_orientation/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' LearnReadOrientationModel"
        " -I {input}"
        " -O {output}"
        " > {log} 2>&1"


rule pileup_summaries:
    input:
        bam = "bqsr/{sample}.bqsr.bam",
        bai = "bqsr/{sample}.bqsr.bai",
    output:
        "contamination/{sample}.pileups.table",
    params:
        java = java_opts,
        common = REF["common_variants"],
    threads: threads_of("pileup_summaries")
    resources:
        mem_mb = mem_of("pileup_summaries"),
    log: "logs/pileup_summaries/{sample}.log"
    benchmark: "benchmarks/pileup_summaries/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' GetPileupSummaries"
        " -I {input.bam}"
        " -V {params.common}"
        " -L {params.common}"
        " -O {output}"
        " > {log} 2>&1"


rule calculate_contamination:
    input:
        "contamination/{sample}.pileups.table",
    output:
        "contamination/{sample}.contamination.table",
    params:
        java = java_opts,
    threads: threads_of("calculate_contamination")
    resources:
        mem_mb = mem_of("calculate_contamination"),
    log: "logs/calculate_contamination/{sample}.log"
    benchmark: "benchmarks/calculate_contamination/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' CalculateContamination"
        " -I {input}"
        " -O {output}"
        " > {log} 2>&1"


rule filter_mutect:
    input:
        vcf = "vcf/raw/{sample}.raw.vcf.gz",
        tbi = "vcf/raw/{sample}.raw.vcf.gz.tbi",
        stats = "vcf/raw/{sample}.raw.vcf.gz.stats",
        rom = "orientation/{sample}.read_orientation_model.tar.gz",
        contamination = "contamination/{sample}.contamination.table",
    output:
        vcf = "vcf/filtered/{sample}.filtered.vcf.gz",
        tbi = "vcf/filtered/{sample}.filtered.vcf.gz.tbi",
    params:
        java = java_opts,
        ref = REF["fasta"],
        extra = TP["filter_mutect"]["extra"],
    threads: threads_of("filter_mutect")
    resources:
        mem_mb = mem_of("filter_mutect"),
    log: "logs/filter_mutect/{sample}.log"
    benchmark: "benchmarks/filter_mutect/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' FilterMutectCalls"
        " -R {params.ref}"
        " -V {input.vcf}"
        " --stats {input.stats}"
        " --ob-priors {input.rom}"
        " --contamination-table {input.contamination}"
        " {params.extra}"
        " -O {output.vcf}"
        " > {log} 2>&1"


rule pon_annotate:
    input:
        vcf = "vcf/filtered/{sample}.filtered.vcf.gz",
        tbi = "vcf/filtered/{sample}.filtered.vcf.gz.tbi",
        pon = PON_VCF if PON_ANNOTATE else [],
    output:
        vcf = "vcf/filtered/{sample}.filtered.pon.vcf.gz",
        tbi = "vcf/filtered/{sample}.filtered.pon.vcf.gz.tbi",
    threads: threads_of("pon_annotate")
    resources:
        mem_mb = mem_of("pon_annotate"),
    log: "logs/pon_annotate/{sample}.log"
    benchmark: "benchmarks/pon_annotate/{sample}.tsv"
    shell:
        "( bcftools annotate -a {input.pon} -m +PON -Oz -o {output.vcf} {input.vcf}"
        " && bcftools index -t {output.vcf} ) > {log} 2>&1"


rule select_pass:
    input:
        vcf = filtered_vcf,
        tbi = lambda wc: f"{filtered_vcf(wc)}.tbi",
    output:
        vcf = "vcf/pass/{sample}.pass.vcf.gz",
        tbi = "vcf/pass/{sample}.pass.vcf.gz.tbi",
    params:
        java = java_opts,
        ref = REF["fasta"],
    threads: threads_of("select_pass")
    resources:
        mem_mb = mem_of("select_pass"),
    log: "logs/select_pass/{sample}.log"
    benchmark: "benchmarks/select_pass/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' SelectVariants"
        " -R {params.ref}"
        " -V {input.vcf}"
        " --exclude-filtered"
        " -O {output.vcf}"
        " > {log} 2>&1"


rule ch_flag:
    input:
        vcf = "vcf/pass/{sample}.pass.vcf.gz",
        tbi = "vcf/pass/{sample}.pass.vcf.gz.tbi",
    output:
        vcf = "vcf/ch/{sample}.ch_flagged.vcf.gz",
        tbi = "vcf/ch/{sample}.ch_flagged.vcf.gz.tbi",
    params:
        java = java_opts,
        ref = REF["fasta"],
        min_dp = config["ch_filter"]["min_dp"],
        min_vaf = config["ch_filter"]["min_vaf"],
    threads: threads_of("ch_flag")
    resources:
        mem_mb = mem_of("ch_flag"),
    log: "logs/ch_flag/{sample}.log"
    benchmark: "benchmarks/ch_flag/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' VariantFiltration"
        " -R {params.ref}"
        " -V {input.vcf}"
        ' --genotype-filter-expression "DP < {params.min_dp}"'
        " --genotype-filter-name LowDP"
        ' --genotype-filter-expression "AF < {params.min_vaf}"'
        " --genotype-filter-name LowVAF"
        " -O {output.vcf}"
        " > {log} 2>&1"


# The genotype filters above only write FORMAT/FT; FILTER stays PASS.
# Filtered genotypes are set to no-call and sites left without a called
# variant genotype are dropped, so the DP/VAF thresholds actually remove calls.
rule ch_select:
    input:
        vcf = "vcf/ch/{sample}.ch_flagged.vcf.gz",
        tbi = "vcf/ch/{sample}.ch_flagged.vcf.gz.tbi",
    output:
        vcf = "vcf/ch/{sample}.ch_final.vcf.gz",
        tbi = "vcf/ch/{sample}.ch_final.vcf.gz.tbi",
    params:
        java = java_opts,
        ref = REF["fasta"],
    threads: threads_of("ch_select")
    resources:
        mem_mb = mem_of("ch_select"),
    log: "logs/ch_select/{sample}.log"
    benchmark: "benchmarks/ch_select/{sample}.tsv"
    shell:
        "gatk --java-options '{params.java}' SelectVariants"
        " -R {params.ref}"
        " -V {input.vcf}"
        " --exclude-filtered"
        " --set-filtered-gt-to-nocall"
        " --exclude-non-variants"
        " -O {output.vcf}"
        " > {log} 2>&1"
