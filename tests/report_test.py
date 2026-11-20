#!/usr/bin/env python3
"""
Test of the CH report on a synthetic cohort, run in CI inside the image:

  python tests/report_test.py [WORKFLOW_HOME]   (default /opt/hemaclone)

Four samples are classified by chip_sample.py and chip_cohort.py, as in a real
run: P1 and P2 with qualifying calls, N1 with a strand-biased call only, and X1
with a qualifying call, listed for exclusion. The report is then rendered with
the R environment of the image (render_report.R). Checked: the HTML and the
four CSV exports exist, X1 is left out of the exports and named in the report,
and the report carries the workflow version.
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
    for folder in ["work", "chip/samples", "logs/batches", "logs/chip_cohort"]:
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

    # ---- settings recorded by the run, as run_batches.py writes them -------------
    settings = {"workflow": {"version": version, "mode": "image", "image": f"hemaclone:{version}"}}
    settings.update(config)
    settings["ch_filter"] = dict(config["ch_filter"], min_dp=100, min_vaf=0.02)
    settings["run"] = {"command": "test", "only": None, "dry_run": False, "snakemake_extra": []}
    with open("logs/batches/run_20000101_000000_settings.yaml", "w") as handle:
        yaml.safe_dump(settings, handle, sort_keys=False)

    with open("work/exclude.txt", "w") as handle:
        handle.write("# left out of the report\n" + "\n".join(EXCLUDED) + "\n")

    # ---- render ------------------------------------------------------------------
    command = [config["tools"]["rscript"], os.path.join(SCRIPTS, "render_report.R"),
               os.path.join(HOME, "workflow", "report", "ch_report.Rmd"),
               "reports/ch_report.html", "reports/results", "work/report_tmp", "work/exclude.txt"]
    os.makedirs("reports", exist_ok=True)
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"FAIL: report not rendered:\n{result.stdout}\n{result.stderr}")

    check(os.path.getsize("reports/ch_report.html") > 0, "reports/ch_report.html written")
    exports = ["ch_samples", "ch_calls", "ch_recurrent_positions", "ch_filter_summary"]
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
    print("Report test passed.")


if __name__ == "__main__":
    main()
