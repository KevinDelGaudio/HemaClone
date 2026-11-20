#!/usr/bin/env python3
"""
Test of the per-sample and cohort CHIP classification on a synthetic sample
with one call per case, run in CI inside the image:

  python tests/chip_sample_test.py [WORKFLOW_HOME]   (default /opt/hemaclone)

Needs bgzip and tabix on PATH. Every new criterion is checked on a call made
to pass or fail it, and the cohort filter counts on the same calls.
"""

import os
import subprocess
import sys
import tempfile
import types

import yaml

HOME = sys.argv[1] if len(sys.argv) > 1 else "/opt/hemaclone"
SCRIPTS = os.path.join(HOME, "workflow", "scripts")

# id, chrom, pos, ref, alt, info flags, SB counts, alt reads, ref reads, gene
CALLS = [
    ("snv_clean",       "chr1", 1000, "C",  "T",  [],      (90, 90, 10, 10), 20, 180, "DNMT3A"),
    ("snv_strand",      "chr1", 1500, "G",  "A",  [],      (85, 85, 30, 0),  30, 170, "TET2"),
    ("indel_giab_low",  "chr1", 2000, "GA", "G",  [],      (97, 97, 3, 3),    6, 194, "ASXL1"),
    ("indel_str_high",  "chr1", 3000, "C",  "CA", ["STR"], (85, 85, 15, 15), 30, 170, "ASXL1"),
    ("indel_plain_low", "chr1", 4000, "TG", "T",  [],      (97, 97, 3, 3),    6, 194, "PPM1D"),
]

# Expected per call: (in_repeat, pass_strand, pass_indel, pass_repeat_indel)
EXPECTED = {
    "snv_clean":       (False, True,  True,  True),
    "snv_strand":      (False, False, True,  True),
    "indel_giab_low":  (True,  True,  False, False),
    "indel_str_high":  (True,  True,  True,  True),
    "indel_plain_low": (False, True,  False, True),
}


def check(condition, message):
    if not condition:
        sys.exit(f"FAIL: {message}")
    print(f"ok    {message}")


def run_script(name, snakemake):
    path = os.path.join(SCRIPTS, name)
    with open(path) as handle:
        code = compile(handle.read(), path, "exec")
    exec(code, {"snakemake": snakemake, "__name__": "__main__"})


def main():
    with open(os.path.join(HOME, "config", "config.yaml")) as handle:
        config = yaml.safe_load(handle)
    chip = dict(config["chip"])
    # Stricter general indel VAF, so that pass_indel is exercised
    chip["indel_min_vaf"] = 0.04

    work = tempfile.mkdtemp()
    os.chdir(work)

    # GIAB-style repeats BED around the indel at chr1:2000 only
    with open("repeats.bed", "w") as handle:
        handle.write("chr1\t1995\t2010\n")
    with open("repeats.bed.gz", "wb") as out:
        subprocess.run(["bgzip", "-c", "repeats.bed"], stdout=out, check=True)
    subprocess.run(["tabix", "-p", "bed", "repeats.bed.gz"], check=True)

    with open("s.vcf", "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n"
                  '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="gnomAD AF">\n'
                  '##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">\n'
                  '##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand counts">\n'
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        for _, chrom, pos, ref, alt, flags, sb, alt_n, ref_n, _ in CALLS:
            info = ";".join(flags) or "."
            vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:AD:SB\t"
                      f"0/1:{ref_n},{alt_n}:{','.join(map(str, sb))}\n")
    subprocess.run(["gzip", "-f", "s.vcf"], check=True)

    with open("s.maf", "w") as maf:
        maf.write("#version 2.4\n")
        maf.write("\t".join(["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
                             "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
                             "Variant_Type", "Protein_Change", "t_alt_count", "t_ref_count",
                             "GNOMAD_AF"]) + "\n")
        for _, chrom, pos, ref, alt, _, _, alt_n, ref_n, gene in CALLS:
            if len(ref) == len(alt):
                vtype, vclass, start = "SNP", "Missense_Mutation", pos
            elif len(ref) > len(alt):
                vtype, vclass, start = "DEL", "Frame_Shift_Del", pos + 1
            else:
                vtype, vclass, start = "INS", "Frame_Shift_Ins", pos
            maf.write("\t".join(map(str, [gene, chrom, start, start, ref, alt, vclass, vtype,
                                          "p.X1Y", alt_n, ref_n, ""])) + "\n")

    sample = types.SimpleNamespace(
        input=types.SimpleNamespace(maf="s.maf", vcf="s.vcf.gz", repeats="repeats.bed.gz",
                                    repeats_tbi="repeats.bed.gz.tbi"),
        output=types.SimpleNamespace(candidates="S1.candidates.tsv", status="S1.status.tsv"),
        log=["S1.log"], wildcards=types.SimpleNamespace(sample="S1"),
        params=types.SimpleNamespace(chip=chip, min_dp=100, min_vaf=0.02, pon_used=True))
    try:
        run_script("chip_sample.py", sample)
    except SystemExit as stop:
        sys.exit(f"FAIL: chip_sample.py stopped: {stop}; log:\n{open('S1.log').read()}")

    import pandas as pd
    table = pd.read_csv("S1.candidates.tsv", sep="\t")
    check(len(table) == len(CALLS), f"{len(CALLS)} calls classified")
    for i, call in enumerate(CALLS):
        name = call[0]
        row = table.iloc[i]
        got = (bool(row["in_repeat"]), bool(row["pass_strand"]),
               bool(row["pass_indel"]), bool(row["pass_repeat_indel"]))
        check(got == EXPECTED[name], f"{name}: in_repeat, pass_strand, pass_indel, pass_repeat_indel = {got}")
    check(abs(table.iloc[0]["sor"] - 0.6931) < 0.0001, f"SOR of balanced counts = ln 2 (got {table.iloc[0]['sor']})")
    check(table.iloc[1]["sor"] > chip["max_sor"], f"SOR of one-strand alt reads above {chip['max_sor']}")
    check(bool(table.iloc[2]["in_giab_repeat"]) and not bool(table.iloc[2]["in_str"]),
          "indel at chr1:2000 found in the repeats BED, not flagged STR")
    check(bool(table.iloc[3]["in_str"]) and not bool(table.iloc[3]["in_giab_repeat"]),
          "indel at chr1:3000 flagged STR, outside the repeats BED")
    check(list(table["qualifies"]) == [True, False, False, True, False], "qualifying calls")

    cohort = types.SimpleNamespace(
        input=types.SimpleNamespace(tables=["S1.candidates.tsv"], status=["S1.status.tsv"]),
        output=types.SimpleNamespace(candidates="c.tsv", maf="c.maf", oncoplot_maf="o.maf",
                                     samples="s.tsv", filter_counts="f.tsv"),
        log=["cohort.log"], params=types.SimpleNamespace(samples=["S1"], large_clone_vaf=0.10))
    run_script("chip_cohort.py", cohort)
    counts = pd.read_csv("f.tsv", sep="\t").set_index("criterion")
    check(counts.loc["pass_strand", "failing"] == 1 and counts.loc["pass_strand", "failing_only_this"] == 1,
          "filter counts: strand fails 1 call, alone")
    check(counts.loc["pass_indel", "failing"] == 2 and counts.loc["pass_indel", "failing_only_this"] == 1,
          "filter counts: indel fails 2 calls, 1 alone")
    check(counts.loc["pass_repeat_indel", "failing"] == 1 and counts.loc["pass_repeat_indel", "failing_only_this"] == 0,
          "filter counts: repeat indel fails 1 call, none alone")
    print("Classification test passed.")


if __name__ == "__main__":
    main()
