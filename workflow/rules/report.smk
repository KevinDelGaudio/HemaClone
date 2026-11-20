# =============================================================================
# Reports, rendered once the cohort tables exist (report block of config.yaml):
# the CH report, the QC report and the methods record. A failed render leaves
# every other output in place; the batch runner lists it, and --cohort-only
# renders it again.
# =============================================================================

rule ch_report:
    input:
        samples = "chip/chip_samples.tsv",
        calls = "chip/chip_calls.maf",
        oncoplot = "chip/chip_calls.oncoplot.maf",
        candidates = "chip/chip_candidates.tsv",
        rmd = REPORT_RMD["ch"],
        # A changed exclusion list renders the report again
        exclude = REPORT["exclude_samples"] or [],
    output:
        html = "reports/ch_report.html",
        samples = "reports/results/ch_samples.csv",
        calls = "reports/results/ch_calls.csv",
        recurrent = "reports/results/ch_recurrent_positions.csv",
        filters = "reports/results/ch_filter_summary.csv",
    params:
        rscript = TOOLS["rscript"],
        script = os.path.join(SCRIPTS, "render_report.R"),
        tmp = os.path.join(config["tmpdir"], "ch_report"),
        exclude = REPORT["exclude_samples"] or "",
    threads: threads_of("ch_report")
    resources:
        mem_mb = mem_of("ch_report"),
    log: "logs/ch_report/ch_report.log"
    benchmark: "benchmarks/ch_report/ch_report.tsv"
    shell:
        "{params.rscript} {params.script} {input.rmd} {output.html} reports/results"
        " {params.tmp} '{params.exclude}' > {log} 2>&1"


rule qc_report:
    input:
        flags = "qc/qc_flags.tsv",
        contamination = expand("contamination/{sample}.contamination.table", sample=SAMPLES),
        relate_samples = "qc/somalier/relate.samples.tsv",
        relate_pairs = "qc/somalier/relate.pairs.tsv",
        rmd = REPORT_RMD["qc"],
        exclude = REPORT["exclude_samples"] or [],
    output:
        html = "reports/qc_report.html",
        summary = "reports/results/qc_summary.csv",
        related = "reports/results/qc_related_pairs.csv",
    params:
        rscript = TOOLS["rscript"],
        script = os.path.join(SCRIPTS, "render_report.R"),
        tmp = os.path.join(config["tmpdir"], "qc_report"),
        exclude = REPORT["exclude_samples"] or "",
    threads: threads_of("qc_report")
    resources:
        mem_mb = mem_of("qc_report"),
    log: "logs/qc_report/qc_report.log"
    benchmark: "benchmarks/qc_report/qc_report.tsv"
    shell:
        "{params.rscript} {params.script} {input.rmd} {output.html} reports/results"
        " {params.tmp} '{params.exclude}' > {log} 2>&1"


# Pinned and installed versions of every package of the image's environments
rule tool_versions:
    input:
        # Re-checked whenever the cohort tables are rebuilt
        samples = "chip/chip_samples.tsv",
    output:
        "reports/results/tool_versions.tsv",
    params:
        docker_dir = os.path.join(WORKFLOW_HOME, "workflow", "docker"),
    threads: threads_of("tool_versions")
    resources:
        mem_mb = mem_of("tool_versions"),
    log: "logs/tool_versions/tool_versions.log"
    benchmark: "benchmarks/tool_versions/tool_versions.tsv"
    script:
        "../scripts/tool_versions.py"


rule methods_report:
    input:
        tools = "reports/results/tool_versions.tsv",
        # Rendered again whenever the cohort tables are rebuilt
        samples = "chip/chip_samples.tsv",
        flags = "qc/qc_flags.tsv",
        rmd = REPORT_RMD["methods"],
    output:
        html = "reports/methods_report.html",
        run = "reports/results/methods_run.csv",
        references = "reports/results/methods_references.csv",
        settings = "reports/results/methods_settings.csv",
    params:
        rscript = TOOLS["rscript"],
        script = os.path.join(SCRIPTS, "render_report.R"),
        tmp = os.path.join(config["tmpdir"], "methods_report"),
        exclude = REPORT["exclude_samples"] or "",
    threads: threads_of("methods_report")
    resources:
        mem_mb = mem_of("methods_report"),
    log: "logs/methods_report/methods_report.log"
    benchmark: "benchmarks/methods_report/methods_report.tsv"
    shell:
        "{params.rscript} {params.script} {input.rmd} {output.html} reports/results"
        " {params.tmp} '{params.exclude}' > {log} 2>&1"
