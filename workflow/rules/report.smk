# =============================================================================
# Reports, rendered once the cohort tables exist (report block of config.yaml).
# A failed render leaves every other output in place; the batch runner lists
# it, and --cohort-only renders it again.
# =============================================================================

rule ch_report:
    input:
        samples = "chip/chip_samples.tsv",
        calls = "chip/chip_calls.maf",
        oncoplot = "chip/chip_calls.oncoplot.maf",
        candidates = "chip/chip_candidates.tsv",
        rmd = CH_RMD,
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
