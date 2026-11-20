#!/usr/bin/env python3
"""
Build the sample sheet (config/samples.tsv) from a folder of paired FASTQ files.

Recognised file name patterns:
  <sample>_S<n>_L<lane>_R<1|2>_001.fastq.gz   lane-split Illumina output
  <sample>_S<n>_R<1|2>_001.fastq.gz           Illumina output without lane splitting
  <sample>_R<1|2>.fastq.gz                    renamed files

Each R1/R2 pair is one "unit". Units sharing the same <sample> are aligned
separately (one read group each) and merged into one BAM at duplicate marking.

BAM and CRAM files in the same folder (<sample>.bam, <sample>.cram), aligned
or not, are one unit each (unit B1): their reads are re-aligned. The sample
name comes from the file name; the library (LB) from the header when it has a
single one. Their rows have the path in the bam column and empty fq1/fq2.

Usage (from the project root):
  python workflow/scripts/make_samplesheet.py --fastq-dir fastq --out config/samples.tsv
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict

PATTERNS = [
    ("lane_split", re.compile(r"^(?P<sample>.+?)_S(?P<snum>\d+)_L(?P<lane>\d+)_R(?P<read>[12])_001\.fastq\.gz$")),
    ("illumina", re.compile(r"^(?P<sample>.+?)_S(?P<snum>\d+)_R(?P<read>[12])_001\.fastq\.gz$")),
    ("plain", re.compile(r"^(?P<sample>.+?)_R(?P<read>[12])\.fastq\.gz$")),
]
NAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")
COLUMNS = ["sample", "unit", "fq1", "fq2", "flowcell", "lane", "barcode", "rg_id", "rg_pu", "lb", "pattern", "bam"]
ALIGNED_SUFFIXES = (".bam", ".cram")
ALIGNED_UNIT = "B1"


def parse_args():
    parser = argparse.ArgumentParser(description="Build the HemaClone sample sheet from a FASTQ folder.")
    parser.add_argument("--fastq-dir", required=True, help="folder with the FASTQ files (relative to the project root)")
    parser.add_argument("--out", required=True, help="output TSV, normally config/samples.tsv")
    parser.add_argument("--force", action="store_true", help="overwrite an existing sample sheet")
    parser.add_argument("--allow-unparsed-headers", action="store_true",
                        help="accept FASTQs whose first header is not Casava 1.8+; read group falls back to <sample>.<unit>")
    return parser.parse_args()


def classify(filename):
    for kind, pattern in PATTERNS:
        match = pattern.match(filename)
        if match:
            return kind, match.groupdict()
    return None, None


def unit_name(kind, fields):
    if kind == "lane_split":
        return f"S{fields['snum']}_L{fields['lane']}"
    if kind == "illumina":
        return f"S{fields['snum']}"
    return "U1"


def first_header(path):
    with gzip.open(path, "rt") as handle:
        line = handle.readline().rstrip("\n")
    if not line.startswith("@"):
        raise ValueError(f"first line is not a FASTQ header: {line[:80]!r}")
    return line[1:]


def read_name(header):
    name = header.split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def parse_casava(header):
    """Return (flowcell, lane, barcode) from a Casava 1.8+ header, or None if the format differs."""
    parts = header.split()
    fields = parts[0].split(":")
    if len(fields) != 7 or not fields[3].isdigit():
        return None
    flowcell = fields[2]
    lane = str(int(fields[3]))
    barcode = "NA"
    if len(parts) > 1:
        comment = parts[1].split(":")
        if len(comment) >= 4 and comment[3] != "":
            barcode = comment[3]
    return flowcell, lane, barcode


def header_read_groups(path):
    """SM, LB and PL values of the @RG lines of a BAM or CRAM header."""
    result = subprocess.run(["samtools", "view", "-H", path], capture_output=True, text=True)
    if result.returncode != 0:
        raise ValueError(result.stderr.strip() or "samtools view -H failed")
    values = {"SM": set(), "LB": set(), "PL": set()}
    for line in result.stdout.splitlines():
        if line.startswith("@RG"):
            for field in line.split("\t")[1:]:
                key, _, value = field.partition(":")
                if key in values:
                    values[key].add(value)
    return values


def aligned_rows(fastq_dir, files, errors, warnings):
    """One sample sheet row per BAM or CRAM file."""
    rows = []
    if files and shutil.which("samtools") is None:
        errors.append("samtools not found on PATH: needed to read BAM/CRAM headers "
                      "(run through hemaclone.sh, which uses the image)")
        return rows
    for filename in files:
        sample, suffix = os.path.splitext(filename)
        path = os.path.join(fastq_dir, filename)
        try:
            groups = header_read_groups(path)
        except ValueError as exc:
            errors.append(f"cannot read the header of {filename}: {exc}")
            continue
        if groups["SM"] and groups["SM"] != {sample}:
            warnings.append(f"{filename}: header SM {sorted(groups['SM'])} differs from the file name; "
                            f"the sample is named {sample}")
        if len(groups["LB"]) == 1:
            library = next(iter(groups["LB"]))
        else:
            library = f"{sample}_lib1"
            if len(groups["LB"]) > 1:
                warnings.append(f"{filename}: {len(groups['LB'])} libraries in the header; treated as one "
                                f"({library}), so duplicates are marked across them")
        rows.append({
            "sample": sample,
            "unit": ALIGNED_UNIT,
            "fq1": "",
            "fq2": "",
            "flowcell": "NA",
            "lane": "NA",
            "barcode": "NA",
            "rg_id": f"{sample}.{ALIGNED_UNIT}",
            "rg_pu": f"{sample}.{ALIGNED_UNIT}",
            "lb": library,
            "pattern": suffix.lstrip("."),
            "bam": path,
        })
    return rows


def main():
    args = parse_args()
    errors = []
    warnings = []

    if not os.path.isdir(args.fastq_dir):
        sys.exit(f"ERROR: FASTQ folder not found: {args.fastq_dir}")
    if os.path.exists(args.out) and not args.force:
        sys.exit(f"ERROR: {args.out} already exists. Review it, or pass --force to overwrite.")

    files = sorted(f for f in os.listdir(args.fastq_dir) if f.endswith((".fastq.gz", ".fq.gz")))
    aligned_files = sorted(f for f in os.listdir(args.fastq_dir) if f.endswith(ALIGNED_SUFFIXES))
    print(f"FASTQ files found: {len(files)}")
    print(f"BAM/CRAM files found: {len(aligned_files)}")

    # (sample, unit) -> {"R1": path, "R2": path, "kind": ..., "file_lane": ...}
    units = defaultdict(dict)
    skipped = []
    for filename in files:
        if filename.startswith("Undetermined"):
            skipped.append(filename)
            continue
        kind, fields = classify(filename)
        if kind is None:
            errors.append(f"unrecognised file name: {filename}")
            continue
        sample = fields["sample"]
        unit = unit_name(kind, fields)
        read = "R" + fields["read"]
        entry = units[(sample, unit)]
        if read in entry:
            errors.append(f"two {read} files for sample {sample}, unit {unit}: {filename}")
            continue
        entry[read] = os.path.join(args.fastq_dir, filename)
        entry["kind"] = kind
        entry["file_lane"] = fields.get("lane")

    if skipped:
        print(f"Skipped Undetermined files: {len(skipped)}")

    # ---- name checks ---------------------------------------------------
    aligned_samples = [os.path.splitext(f)[0] for f in aligned_files]
    for sample in sorted({s for s, _ in units} & set(aligned_samples)):
        errors.append(f"sample {sample} has both FASTQ and BAM/CRAM files")
    for sample in sorted({s for s in aligned_samples if aligned_samples.count(s) > 1}):
        errors.append(f"sample {sample} has both a BAM and a CRAM file")
    sample_names = sorted({sample for sample, _ in units} | set(aligned_samples))
    for sample in sample_names:
        if not NAME_PATTERN.match(sample):
            errors.append(f"sample name has characters outside [A-Za-z0-9_-]: {sample}")
    lowered = defaultdict(list)
    for sample in sample_names:
        lowered[sample.lower()].append(sample)
    for names in lowered.values():
        if len(names) > 1:
            errors.append(f"sample names differ only by case: {names}")

    # ---- pairs and headers ---------------------------------------------
    rows = []
    for (sample, unit), entry in sorted(units.items()):
        if "R1" not in entry or "R2" not in entry:
            present = entry.get("R1") or entry.get("R2")
            errors.append(f"missing mate for sample {sample}, unit {unit}: only {present}")
            continue

        try:
            header1 = first_header(entry["R1"])
            header2 = first_header(entry["R2"])
        except (OSError, ValueError) as exc:
            errors.append(f"cannot read header of sample {sample}, unit {unit}: {exc}")
            continue

        if read_name(header1) != read_name(header2):
            errors.append(f"R1/R2 first read names differ for sample {sample}, unit {unit}: "
                          f"{read_name(header1)} vs {read_name(header2)}")
            continue

        parsed = parse_casava(header1)
        if parsed is None:
            if not args.allow_unparsed_headers:
                errors.append(f"header not in Casava 1.8+ format for sample {sample}, unit {unit}: {header1[:80]}")
                continue
            warnings.append(f"unparsed header, read group falls back to {sample}.{unit}: {header1[:80]}")
            flowcell, lane, barcode = "NA", "NA", "NA"
            rg_id = f"{sample}.{unit}"
        else:
            flowcell, lane, barcode = parsed
            if entry["file_lane"] is not None and int(entry["file_lane"]) != int(lane):
                errors.append(f"file name lane L{entry['file_lane']} differs from header lane {lane} "
                              f"for sample {sample}, unit {unit}")
                continue
            rg_id = f"{flowcell}.{lane}"

        rows.append({
            "sample": sample,
            "unit": unit,
            "fq1": entry["R1"],
            "fq2": entry["R2"],
            "flowcell": flowcell,
            "lane": lane,
            "barcode": barcode,
            "rg_id": rg_id,
            "rg_pu": f"{flowcell}.{lane}.{barcode}",
            "lb": f"{sample}_lib1",
            "pattern": entry["kind"],
            "bam": "",
        })

    rows += aligned_rows(args.fastq_dir, aligned_files, errors, warnings)

    # ---- read group IDs must be unique within a sample ----------------
    seen_rg = defaultdict(list)
    for row in rows:
        seen_rg[(row["sample"], row["rg_id"])].append(row["unit"])
    for (sample, rg_id), unit_list in seen_rg.items():
        if len(unit_list) > 1:
            errors.append(f"read group ID {rg_id} shared by units {unit_list} of sample {sample}")

    if errors:
        print(f"\n{len(errors)} ERROR(S), sample sheet NOT written:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # ---- summary -------------------------------------------------------
    units_per_sample = defaultdict(list)
    pattern_counts = defaultdict(int)
    for row in rows:
        units_per_sample[row["sample"]].append(row["unit"])
        pattern_counts[row["pattern"]] += 1

    print(f"\nSample sheet written: {args.out}")
    print(f"  samples : {len(units_per_sample)}")
    print(f"  units   : {len(rows)}")
    for kind, count in sorted(pattern_counts.items()):
        print(f"  units with {kind} names: {count}")
    multi = {s: u for s, u in units_per_sample.items() if len(u) > 1}
    print(f"\nSamples with more than one unit (merged at duplicate marking): {len(multi)}")
    for sample, unit_list in sorted(multi.items()):
        print(f"  {sample}: {', '.join(unit_list)}")
    if warnings:
        print(f"\n{len(warnings)} warning(s):")
        for warning in warnings:
            print(f"  - {warning}")
    print(f"\nTo review it: column -t {args.out} | less -S")


if __name__ == "__main__":
    main()
