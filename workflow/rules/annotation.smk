# =============================================================================
# gnomAD allele frequency, then Funcotator, two tiers per sample:
#   maf/{sample}_all.funcotated.maf       PASS calls, before CH DP/VAF filters
#   maf/{sample}_ch_final.funcotated.maf  CH call set after DP/VAF filters
#
# INFO/GNOMAD_AF is copied from the local gnomAD VCF (reference.gnomad, the
# same file Mutect2 uses as germline resource), matched by position and
# allele. Funcotator carries INFO fields into the MAF, so every MAF has a
# GNOMAD_AF column; a call absent from gnomAD has it empty.
# =============================================================================

TIER_VCF = {
    "all": "vcf/pass/{sample}.pass.vcf.gz",
    "ch_final": "vcf/ch/{sample}.ch_final.vcf.gz",
}


rule gnomad_annotate:
    input:
        vcf = lambda wc: TIER_VCF[wc.tier].format(sample=wc.sample),
        tbi = lambda wc: f"{TIER_VCF[wc.tier].format(sample=wc.sample)}.tbi",
        gnomad = REF["gnomad"],
    output:
        vcf = "vcf/gnomad/{sample}.{tier}.vcf.gz",
        tbi = "vcf/gnomad/{sample}.{tier}.vcf.gz.tbi",
    threads: threads_of("gnomad_annotate")
    resources:
        mem_mb = mem_of("gnomad_annotate"),
    log: "logs/gnomad_annotate/{sample}_{tier}.log"
    benchmark: "benchmarks/gnomad_annotate/{sample}_{tier}.tsv"
    shell:
        "( bcftools annotate -a {input.gnomad} -c INFO/GNOMAD_AF:=INFO/AF"
        " -Oz -o {output.vcf} {input.vcf}"
        " && bcftools index -t {output.vcf} ) > {log} 2>&1"


rule funcotator:
    input:
        vcf = "vcf/gnomad/{sample}.{tier}.vcf.gz",
        tbi = "vcf/gnomad/{sample}.{tier}.vcf.gz.tbi",
    output:
        "maf/{sample}_{tier}.funcotated.maf",
    params:
        java = java_opts,
        ref = REF["fasta"],
        sources = REF["funcotator_ds"],
        ref_version = TP["funcotator"]["ref_version"],
        extra = TP["funcotator"]["extra"],
    threads: threads_of("funcotator")
    resources:
        mem_mb = mem_of("funcotator"),
    log: "logs/funcotator/{sample}_{tier}.log"
    benchmark: "benchmarks/funcotator/{sample}_{tier}.tsv"
    shell:
        "gatk --java-options '{params.java}' Funcotator"
        " --variant {input.vcf}"
        " --reference {params.ref}"
        " --ref-version {params.ref_version}"
        " --data-sources-path {params.sources}"
        " --output {output}"
        " --output-file-format MAF"
        " {params.extra}"
        " > {log} 2>&1"
