# Per-sample PASS/FLAG table from Picard HsMetrics, InsertSize and
# MarkDuplicates metrics. Called by the qc_flags rule (Snakemake script).

import csv
import sys

COLUMNS = ["sample", "mean_target_coverage", "pct_target_20x", "pct_target_50x",
           "pct_selected_bases", "fold_enrichment", "pct_exc_dupe", "at_dropout",
           "gc_dropout", "median_insert_size", "pct_dup_markdup", "status", "flags"]


def picard_metrics(path):
    """First metrics row of a Picard metrics file as a dict; empty dict if absent."""
    try:
        with open(path) as handle:
            lines = handle.read().splitlines()
    except FileNotFoundError:
        return {}
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            if i + 2 >= len(lines):
                return {}
            header = lines[i + 1].split("\t")
            values = lines[i + 2].split("\t")
            return dict(zip(header, values))
    return {}


def value(metrics, column):
    text = metrics.get(column, "")
    return text if text != "" else "NA"


def as_float(text):
    try:
        return float(text)
    except ValueError:
        return None


def flags_for(row, thr):
    flags = []
    cov = as_float(row["mean_target_coverage"])
    t20 = as_float(row["pct_target_20x"])
    sel = as_float(row["pct_selected_bases"])
    dupe = as_float(row["pct_exc_dupe"])
    atd = as_float(row["at_dropout"])
    gcd = as_float(row["gc_dropout"])

    if cov is not None and cov < thr["min_mean_coverage"]:
        flags.append("LOW_COVERAGE")
    if t20 is not None and t20 < thr["min_pct_target_20x"]:
        flags.append("LOW_20X_BREADTH")
    if sel is not None and sel < thr["min_pct_selected"]:
        flags.append("POOR_CAPTURE")
    if dupe is not None and dupe > thr["max_pct_dupe"]:
        flags.append("HIGH_DUPLICATION")
    if atd is not None and atd > thr["max_dropout"]:
        flags.append("AT_DROPOUT")
    if gcd is not None and gcd > thr["max_dropout"]:
        flags.append("GC_DROPOUT")
    if cov is None:
        flags.append("NO_HS_METRICS")
    return flags


def main(snakemake):
    log = open(snakemake.log[0], "w")
    samples = snakemake.params.samples
    thr = snakemake.params.thresholds
    hs_files = list(snakemake.input.hs)
    insert_files = list(snakemake.input.insert_size)
    dup_files = list(snakemake.input.dup)

    if not (len(samples) == len(hs_files) == len(insert_files) == len(dup_files)):
        sys.exit("qc_flags: sample and metrics file counts differ")

    rows = []
    for sample, hs_path, insert_path, dup_path in zip(samples, hs_files, insert_files, dup_files):
        hs = picard_metrics(hs_path)
        ins = picard_metrics(insert_path)
        dup = picard_metrics(dup_path)
        row = {
            "sample": sample,
            "mean_target_coverage": value(hs, "MEAN_TARGET_COVERAGE"),
            "pct_target_20x": value(hs, "PCT_TARGET_BASES_20X"),
            "pct_target_50x": value(hs, "PCT_TARGET_BASES_50X"),
            "pct_selected_bases": value(hs, "PCT_SELECTED_BASES"),
            "fold_enrichment": value(hs, "FOLD_ENRICHMENT"),
            "pct_exc_dupe": value(hs, "PCT_EXC_DUPE"),
            "at_dropout": value(hs, "AT_DROPOUT"),
            "gc_dropout": value(hs, "GC_DROPOUT"),
            "median_insert_size": value(ins, "MEDIAN_INSERT_SIZE"),
            "pct_dup_markdup": value(dup, "PERCENT_DUPLICATION"),
        }
        flags = flags_for(row, thr)
        row["status"] = "FLAG" if flags else "PASS"
        row["flags"] = ";".join(flags) if flags else "-"
        rows.append(row)

    with open(snakemake.output[0], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    flagged = [r for r in rows if r["status"] == "FLAG"]
    print(f"Samples evaluated: {len(rows)}", file=log)
    print(f"Flagged: {len(flagged)}", file=log)
    for r in flagged:
        print(f"  {r['sample']}: cov={r['mean_target_coverage']} 20x={r['pct_target_20x']} "
              f"dup={r['pct_exc_dupe']} {r['flags']}", file=log)
    log.close()


main(snakemake)
