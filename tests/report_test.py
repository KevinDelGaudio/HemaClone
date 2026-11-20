#!/usr/bin/env python3
"""
Test of the three reports on a synthetic cohort, run in CI inside the image:

  python tests/report_test.py [WORKFLOW_HOME]   (default /opt/hemaclone)

Four samples are classified by chip_sample.py and chip_cohort.py, as in a real
run: P1 and P2 with qualifying calls, N1 with a strand-biased call only, and X1
with a qualifying call, listed for exclusion. QC tables (qc_flags.tsv,
contamination, somalier pairs), a sample sheet and a resources manifest are
written as the workflow writes them, and tool_versions.py runs as in a run.
The reports are then rendered with the R environment of the image
(render_report.R). Checked: every HTML and CSV export exists; X1 is left out of
the CH exports, named in the CH report and marked excluded in the QC summary;
the reports carry the workflow version; the methods record lists the settings
and finds the reference files in the manifest.
"""

import csv
import os
import subprocess
import sys
import tempfile
import types

import yaml

HOME = sys.argv[1] if len(sys.argv) > 1 else "/opt/hemaclone"
SCRIPTS = os.path.join(HOME, "workflow", "scripts")

# sample -> calls: (chrom, pos, ref, alt, info flags, SB counts, alt reads, ref reads,
#                   gene, Variant_Classification, Protein_Change)
COHORT = {
    "P1": [("chr1", 1000, "C", "T", [], (90, 90, 10, 10), 20, 180, "DNMT3A", "Missense_Mutation", "p.R882H"),
           ("chr1", 1500, "G", "A", [], (85, 85, 15, 15), 30, 170, "TET2", "Nonsense_Mutation", "p.Q1034*")],
    "P2": [("chr1", 1000, "G", "A", [], (92, 93, 8, 7), 15, 185, "DNMT3A", "Missense_Mutation", "p.R882C"),
           ("chr1", 3000, "C", "CA", ["STR"], (85, 85, 15, 15), 30, 170, "ASXL1", "Frame_Shift_Ins", "p.G646Wfs*12")],
    "N1": [("chr1", 1500, "G", "A", [], (85, 85, 30, 0), 30, 170, "TET2", "Nonsense_Mutation", "p.Q1034*")],
    "X1": [("chr1", 1000, "C", "T", [], (90, 90, 10, 10), 20, 180, "DNMT3A", "Missense_Mutation", "p.R882H")],
}
EXCLUDED = ["X1"]


def check(condition, message):
    if not condition:
        sys.exit(f"FAIL: {message}")
    print(f"ok    {message}")


def run_script(name, snakemake):
    path = os.path.join(SCRIPTS, name)
    with open(path) as handle:
        code = compile(handle.read(), path, "exec")
    exec(code, {"snakemake": snakemake, "__name__": "__main__"})


def write_sample_inputs(sample, calls):
    """Mutect2-style VCF and Funcotator-style MAF of one sample, in work/."""
    with open(f"work/{sample}.vcf", "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n"
                  '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD AF">\n'
                  '##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">\n'
                  '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand counts">\n'
                  f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for chrom, pos, ref, alt, flags, sb, alt_n, ref_n, _, _, _ in calls:
            info = ";".join(flags) or "."
            vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:AD:SB\t"
                      f"0/1:{ref_n},{alt_n}:{','.join(map(str, sb))}\n")
    subprocess.run(["gzip", "-f", f"work/{sample}.vcf"], check=True)

    with open(f"work/{sample}.maf", "w") as maf:
        maf.write("#version 2.4\n")
        maf.write("\t".join(["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                             "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
                             "Variant_Type", "Protein_Change", "t_alt_count", "t_ref_count",
                             "GNOMAD_AF"]) + "\n")
        for chrom, pos, ref, alt, _, _, alt_n, ref_n, gene, vclass, protein in calls:
            if len(ref) == len(alt):
                vtype, start = "SNP", pos
            elif len(ref) > len(alt):
                vtype, start = "DEL", pos + 1
            else:
                vtype, start = "INS", pos
            maf.write("\t".join(map(str, [gene, chrom, start, start, ref, alt, vclass, vtype,
                                          protein, alt_n, ref_n, ""])) + "\n")


def read_csv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def main():
    with open(os.path.join(HOME, "config", "config.yaml")) as handle:
        config = yaml.safe_load(handle)
    with open(os.path.join(HOME, "VERSION")) as handle:
        version = handle.read().strip()

    project = tempfile.mkdtemp()
    os.chdir(project)
    for folder in ["work", "chip/samples", "logs/batches", "logs/chip_cohort", "qc/somalier",
                   "contamination", "config", "resources", "reports/results"]:
        os.makedirs(folder, exist_ok=True)

    # Empty repeats BED: no call is in a GIAB repeat
    open("work/repeats.bed", "w").close()
    with open("work/repeats.bed.gz", "wb") as out:
        subprocess.run(["bgzip", "-c", "work/repeats.bed"], stdout=out, check=True)
    subprocess.run(["tabix", "-p", "bed", "work/repeats.bed.gz"], check=True)

    # ---- per-sample and cohort classification, as in a run ------------------------
    samples = list(COHORT)
    for sample, calls in COHORT.items():
        write_sample_inputs(sample, calls)
        job = types.SimpleNamespace(
            input=types.SimpleNamespace(maf=f"work/{sample}.maf", vcf=f"work/{sample}.vcf.gz",
                                        repeats="work/repeats.bed.gz",
                                        repeats_tbi="work/repeats.bed.gz.tbi"),
            output=types.SimpleNamespace(candidates=f"chip/samples/{sample}.candidates.tsv",
                                         status=f"chip/samples/{sample}.status.tsv"),
            log=[f"work/{sample}.log"], wildcards=types.SimpleNamespace(sample=sample),
            params=types.SimpleNamespace(chip=config["chip"], min_dp=100, min_vaf=0.02, pon_used=True))
        run_script("chip_sample.py", job)

    cohort = types.SimpleNamespace(
        input=types.SimpleNamespace(tables=[f"chip/samples/{s}.candidates.tsv" for s in samples],
                                    status=[f"chip/samples/{s}.status.tsv" for s in samples]),
        output=types.SimpleNamespace(candidates="chip/chip_candidates.tsv", maf="chip/chip_calls.maf",
                                     oncoplot_maf="chip/chip_calls.oncoplot.maf",
                                     samples="chip/chip_samples.tsv",
                                     filter_counts="chip/chip_filter_counts.tsv"),
        log=["logs/chip_cohort/chip_cohort.log"],
        params=types.SimpleNamespace(samples=samples, large_clone_vaf=config["chip"]["large_clone_vaf"]))
    run_script("chip_cohort.py", cohort)

    status = {row["Tumor_Sample_Barcode"]: row["CHIP_status"]
              for row in csv.DictReader(open("chip/chip_samples.tsv"), delimiter="\t")}
    check(status == {"P1": "CHIP+", "P2": "CHIP+", "N1": "CHIP-", "X1": "CHIP+"},
          f"cohort classified: {status}")

    # ---- QC outputs, sample sheet and resources manifest, as the workflow writes them
    sheet_columns = ["sample", "unit", "fq1", "fq2", "flowcell", "lane", "barcode",
                     "rg_id", "rg_pu", "lb", "pattern", "bam"]
    with open("config/samples.tsv", "w") as handle:
        handle.write("\t".join(sheet_columns) + "\n")
        for sample in samples:
            if sample == "N1":
                row = [sample, "B1", "", "", "NA", "NA", "NA", f"{sample}.B1", f"{sample}.B1",
                       f"{sample}_lib1", "cram", f"fastq/{sample}.cram"]
            else:
                row = [sample, "U1", f"fastq/{sample}_R1.fastq.gz", f"fastq/{sample}_R2.fastq.gz",
                       "FC1", "1", "NA", "FC1.1", "FC1.1.NA", f"{sample}_lib1", "plain", ""]
            handle.write("\t".join(row) + "\n")

    flag_columns = ["sample", "mean_target_coverage", "pct_target_20x", "pct_target_50x",
                    "pct_selected_bases", "fold_enrichment", "pct_exc_dupe", "at_dropout",
                    "gc_dropout", "median_insert_size", "pct_dup_markdup", "x_depth_ratio",
                    "y_depth_ratio", "inferred_sex", "related_samples", "status", "flags"]
    qc_values = {
        "P1": ["150", "0.97", "0.90", "0.70", "40", "0.10", "1.2", "0.8", "210", "0.12", "1.01", "0.002",
               "NA", "-", "PASS", "-"],
        "P2": ["140", "0.96", "0.88", "0.71", "41", "0.11", "1.1", "0.9", "205", "0.13", "0.52", "0.48",
               "NA", "X1:0.98", "FLAG", "POSSIBLE_DUPLICATE"],
        "N1": ["25", "0.60", "0.20", "0.69", "39", "0.55", "1.5", "1.2", "180", "0.58", "1.00", "0.001",
               "NA", "-", "FLAG", "LOW_COVERAGE;LOW_20X_BREADTH;HIGH_DUPLICATION"],
        "X1": ["138", "0.95", "0.87", "0.70", "40", "0.11", "1.1", "0.9", "206", "0.12", "0.51", "0.47",
               "NA", "P2:0.98", "FLAG", "POSSIBLE_DUPLICATE"],
    }
    with open("qc/qc_flags.tsv", "w") as handle:
        handle.write("\t".join(flag_columns) + "\n")
        for sample in samples:
            handle.write("\t".join([sample] + qc_values[sample]) + "\n")
    for sample, fraction in zip(samples, ["0.001", "0.002", "0.031", "0.002"]):
        with open(f"contamination/{sample}.contamination.table", "w") as handle:
            handle.write(f"sample\tcontamination\terror\n{sample}\t{fraction}\t0.0005\n")
    with open("qc/somalier/relate.pairs.tsv", "w") as handle:
        handle.write("#sample_a\tsample_b\trelatedness\tibs0\tibs2\tn\n")
        for a, b, value in [("P1", "P2", "0.01"), ("P1", "N1", "-0.02"), ("P1", "X1", "0.02"),
                            ("P2", "N1", "0.00"), ("P2", "X1", "0.98"), ("N1", "X1", "0.01")]:
            handle.write(f"{a}\t{b}\t{value}\t10\t500\t2000\n")

    resources = os.path.join(project, "resources")
    with open(os.path.join(resources, "MANIFEST.tsv"), "w") as handle:
        handle.write("file\tsource\tbytes\tmd5\tdate\n")
        handle.write(f"{config['reference']['fasta']}\thttps://example.org/{config['reference']['fasta']}"
                     "\t100\t0123456789abcdef0123456789abcdef\t2026-01-01 00:00\n")

    tools = types.SimpleNamespace(
        output=["reports/results/tool_versions.tsv"], log=["work/tool_versions.log"],
        params=types.SimpleNamespace(docker_dir=os.path.join(HOME, "workflow", "docker")))
    run_script("tool_versions.py", tools)
    tool_rows = list(csv.DictReader(open("reports/results/tool_versions.tsv"), delimiter="\t"))
    check(any(r["package"] == "bioconductor-maftools" for r in tool_rows),
          f"tool_versions.tsv lists the pinned packages ({len(tool_rows)})")

    # ---- settings recorded by the run, as run_batches.py writes them -------------
    settings = {"workflow": {"version": version, "mode": "image", "image": f"hemaclone:{version}",
                             "code": HOME}}
    settings.update(config)
    settings["ch_filter"] = dict(config["ch_filter"], min_dp=100, min_vaf=0.02)
    settings["reference"] = dict(config["reference"], dir=resources)
    settings["samples"] = "config/samples.tsv"
    settings["machine"] = {"cores": 4, "mem_gb": 16, "schedule_mb": 14745}
    settings["reference_resolved"] = {k: (v if k == "dir" else os.path.join(resources, v))
                                      for k, v in settings["reference"].items()}
    settings["pon_resolved"] = os.path.join(resources, config["pon"]["vcf"])
    settings["run"] = {"command": "test", "only": None, "dry_run": False, "snakemake_extra": []}
    with open("logs/batches/run_20000101_000000_settings.yaml", "w") as handle:
        yaml.safe_dump(settings, handle, sort_keys=False)

    with open("work/exclude.txt", "w") as handle:
        handle.write("# left out of the report\n" + "\n".join(EXCLUDED) + "\n")

    # ---- render ------------------------------------------------------------------
    for name in ["ch_report", "qc_report", "methods_report"]:
        command = [config["tools"]["rscript"], os.path.join(SCRIPTS, "render_report.R"),
                   os.path.join(HOME, "workflow", "report", f"{name}.Rmd"),
                   f"reports/{name}.html", "reports/results", f"work/{name}_tmp", "work/exclude.txt"]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(f"FAIL: {name} not rendered:\n{result.stdout}\n{result.stderr}")
        check(os.path.getsize(f"reports/{name}.html") > 0, f"reports/{name}.html written")

    exports = ["ch_samples", "ch_calls", "ch_recurrent_positions", "ch_filter_summary",
               "qc_summary", "qc_related_pairs", "methods_run", "methods_references", "methods_settings"]
    for name in exports:
        check(os.path.isfile(f"reports/results/{name}.csv"), f"reports/results/{name}.csv written")

    rows = read_csv("reports/results/ch_samples.csv")
    ids = sorted(row["Tumor_Sample_Barcode"] for row in rows)
    check(ids == ["N1", "P1", "P2"], f"excluded sample left out of ch_samples.csv: {ids}")
    calls = read_csv("reports/results/ch_calls.csv")
    check(len(calls) == 4 and all(row["Tumor_Sample_Barcode"] != "X1" for row in calls),
          f"4 qualifying calls exported, none from X1 (got {len(calls)})")

    html = open("reports/ch_report.html").read()
    check("Samples left out of this report (1):" in html and "X1" in html, "excluded sample named in the report")
    check(f"HemaClone {version}" in html, f"report carries the workflow version {version}")
    check(os.path.isfile(os.path.join(project, "chip", "chip_samples.tsv")) and
          "X1" in open("chip/chip_samples.tsv").read(), "chip/ tables keep the excluded sample")

    summary = {row["sample"]: row for row in read_csv("reports/results/qc_summary.csv")}
    check(sorted(summary) == sorted(samples), "QC summary keeps every sample")
    check(summary["X1"]["excluded"] == "TRUE" and summary["P1"]["excluded"] == "FALSE",
          "QC summary marks the excluded sample")
    check(summary["N1"]["input"] == "CRAM" and summary["P1"]["input"] == "FASTQ", "QC summary: input type")
    pairs = read_csv("reports/results/qc_related_pairs.csv")
    check(len(pairs) == 1 and {pairs[0]["sample_a"], pairs[0]["sample_b"]} == {"P2", "X1"},
          "related pairs: P2 and X1 only")
    for name in ["qc_report", "methods_report"]:
        check(f"HemaClone {version}" in open(f"reports/{name}.html").read(), f"{name} carries the version")
    settings_rows = {row["setting"]: row["value"] for row in read_csv("reports/results/methods_settings.csv")}
    check(settings_rows.get("chip.max_sor") == str(config["chip"]["max_sor"]),
          "methods record lists chip.max_sor")
    references = {row["setting"]: row["in_manifest"] for row in read_csv("reports/results/methods_references.csv")}
    check(references.get("reference.fasta") == "TRUE" and references.get("reference.dbsnp") == "FALSE",
          "methods record matches reference files to the manifest")
    print("Report test passed.")


if __name__ == "__main__":
    main()
