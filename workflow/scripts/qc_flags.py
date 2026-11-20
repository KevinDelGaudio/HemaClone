# Per-sample PASS/FLAG table from Picard coverage (HsMetrics for WES,
# WgsMetrics for WGS), InsertSize and MarkDuplicates metrics, and the somalier
# sex and identity check. Called by the qc_flags rule (Snakemake script).

import csv
import sys

IDENTITY_COLUMNS = ["x_depth_ratio", "y_depth_ratio", "inferred_sex", "related_samples"]

COLUMNS = ["sample", "mean_target_coverage", "pct_target_20x", "pct_target_50x",
           "pct_selected_bases", "fold_enrichment", "pct_exc_dupe", "at_dropout",
           "gc_dropout", "median_insert_size", "pct_dup_markdup"] + IDENTITY_COLUMNS + ["status", "flags"]

WGS_COLUMNS = ["sample", "mean_coverage", "pct_20x", "pct_exc_dupe", "pct_exc_mapq",
               "median_insert_size", "pct_dup_markdup"] + IDENTITY_COLUMNS + ["status", "flags"]


def fail(log, message):
    """Write the error to the rule log (Snakemake does not capture a script's
    stderr there), then stop."""
    print(f"ERROR: {message}", file=log, flush=True)
    sys.exit(message)


def read_tsv(path, log):
    """Header (leading # removed) and rows of a somalier TSV."""
    with open(path) as handle:
        lines = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not lines:
        fail(log, f"{path}: empty")
    header = [name.lstrip("#") for name in lines[0]]
    return header, [dict(zip(header, row)) for row in lines[1:]]


def somalier_ratios(path, samples, log):
    """Sample -> (X/autosome, Y/autosome) mean depth at the somalier sites."""
    header, rows = read_tsv(path, log)
    needed = ["sample_id", "gt_depth_mean", "X_depth_mean", "Y_depth_mean"]
    missing = [c for c in needed if c not in header]
    if missing:
        fail(log, f"{path}: columns {missing} not found; header: {header}")
    ratios = {}
    for row in rows:
        autosomal = as_float(row["gt_depth_mean"])
        x, y = as_float(row["X_depth_mean"]), as_float(row["Y_depth_mean"])
        if autosomal and x is not None and y is not None:
            ratios[row["sample_id"]] = (round(x / autosomal, 4), round(y / autosomal, 4))
        else:
            ratios[row["sample_id"]] = (None, None)
    absent = [s for s in samples if s not in ratios]
    if absent:
        fail(log, f"{path}: samples missing from the somalier output: {absent}")
    return ratios


def related_pairs(path, minimum, log):
    """Sample -> list of (other sample, relatedness) with relatedness >= minimum."""
    header, rows = read_tsv(path, log)
    needed = ["sample_a", "sample_b", "relatedness"]
    missing = [c for c in needed if c not in header]
    if missing:
        fail(log, f"{path}: columns {missing} not found; header: {header}")
    related = {}
    for row in rows:
        value = as_float(row["relatedness"])
        if value is None or value < minimum:
            continue
        related.setdefault(row["sample_a"], []).append((row["sample_b"], value))
        related.setdefault(row["sample_b"], []).append((row["sample_a"], value))
    return related


def read_expected_sex(path, log):
    expected = {}
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.strip().split("\t")
            if not line.strip() or line.startswith("#") or fields[0] == "sample":
                continue
            if len(fields) < 2 or fields[1].lower() not in ("male", "female"):
                fail(log, f"{path} line {number}: expected 'sample<TAB>male|female'")
            expected[fields[0]] = fields[1].lower()
    return expected


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


def flags_for_wgs(row, thr):
    flags = []
    cov = as_float(row["mean_coverage"])
    t20 = as_float(row["pct_20x"])
    dupe = as_float(row["pct_exc_dupe"])
    if cov is not None and cov < thr["min_mean_coverage"]:
        flags.append("LOW_COVERAGE")
    if t20 is not None and t20 < thr["min_pct_20x"]:
        flags.append("LOW_20X_BREADTH")
    if dupe is not None and dupe > thr["max_pct_dupe"]:
        flags.append("HIGH_DUPLICATION")
    if cov is None:
        flags.append("NO_WGS_METRICS")
    return flags


def wgs_row(sample, coverage, ins, dup):
    return {
        "sample": sample,
        "mean_coverage": value(coverage, "MEAN_COVERAGE"),
        "pct_20x": value(coverage, "PCT_20X"),
        "pct_exc_dupe": value(coverage, "PCT_EXC_DUPE"),
        "pct_exc_mapq": value(coverage, "PCT_EXC_MAPQ"),
        "median_insert_size": value(ins, "MEDIAN_INSERT_SIZE"),
        "pct_dup_markdup": value(dup, "PERCENT_DUPLICATION"),
    }


def main(snakemake, log):
    samples = snakemake.params.samples
    thr = snakemake.params.thresholds
    wgs = snakemake.params.data_type == "wgs"
    hs_files = list(snakemake.input.coverage)
    insert_files = list(snakemake.input.insert_size)
    dup_files = list(snakemake.input.dup)

    if not (len(samples) == len(hs_files) == len(insert_files) == len(dup_files)):
        fail(log, "qc_flags: sample and metrics file counts differ")

    identity = snakemake.params.identity
    ratios = somalier_ratios(snakemake.input.relate_samples, samples, log)
    related = related_pairs(snakemake.input.relate_pairs, identity["min_relatedness_list"], log)
    expected = read_expected_sex(identity["expected_sex"], log) if identity["expected_sex"] else {}
    if expected:
        not_listed = [s for s in samples if s not in expected]
        if not_listed:
            print(f"Expected sex not given for {len(not_listed)} sample(s): {', '.join(not_listed)}", file=log)

    def add_identity(row, flags):
        x_ratio, y_ratio = ratios[row["sample"]]
        row["x_depth_ratio"] = "NA" if x_ratio is None else x_ratio
        row["y_depth_ratio"] = "NA" if y_ratio is None else y_ratio
        sex = "NA"
        if identity["sex_y_ratio_male"] is not None and y_ratio is not None:
            sex = "male" if y_ratio >= identity["sex_y_ratio_male"] else "female"
        row["inferred_sex"] = sex
        pairs = sorted(related.get(row["sample"], []), key=lambda p: -p[1])
        row["related_samples"] = ",".join(f"{other}:{value:.2f}" for other, value in pairs) or "-"
        if any(value >= identity["min_relatedness_flag"] for _, value in pairs):
            flags.append("POSSIBLE_DUPLICATE")
        if sex != "NA" and row["sample"] in expected and expected[row["sample"]] != sex:
            flags.append("SEX_MISMATCH")

    rows = []
    for sample, hs_path, insert_path, dup_path in zip(samples, hs_files, insert_files, dup_files):
        hs = picard_metrics(hs_path)
        ins = picard_metrics(insert_path)
        dup = picard_metrics(dup_path)
        if wgs:
            row = wgs_row(sample, hs, ins, dup)
            flags = flags_for_wgs(row, thr)
            add_identity(row, flags)
            row["status"] = "FLAG" if flags else "PASS"
            row["flags"] = ";".join(flags) if flags else "-"
            rows.append(row)
            continue
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
        add_identity(row, flags)
        row["status"] = "FLAG" if flags else "PASS"
        row["flags"] = ";".join(flags) if flags else "-"
        rows.append(row)

    with open(snakemake.output[0], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=WGS_COLUMNS if wgs else COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    flagged = [r for r in rows if r["status"] == "FLAG"]
    print(f"Samples evaluated: {len(rows)}", file=log)
    print(f"Flagged: {len(flagged)}", file=log)
    coverage_key, breadth_key = ("mean_coverage", "pct_20x") if wgs else ("mean_target_coverage", "pct_target_20x")
    for r in flagged:
        print(f"  {r['sample']}: cov={r[coverage_key]} 20x={r[breadth_key]} "
              f"dup={r['pct_exc_dupe']} {r['flags']}", file=log)


with open(snakemake.log[0], "w") as log_handle:
    main(snakemake, log_handle)
