# Cohort CHIP tables from the per-sample classifications (Snakemake script).
#
# Outputs:
#   chip/chip_candidates.tsv    every ch_final call of every sample, all criteria
#   chip/chip_calls.maf         qualifying calls counted for CHIP status (maftools)
#   chip/chip_calls.oncoplot.maf  the same plus one Silent placeholder row per
#                               CHIP- sample, so oncoplot shows the whole cohort
#   chip/chip_samples.tsv       one row per sample, CHIP- included (clinicalData)
#   chip/chip_filter_counts.tsv per criterion: calls failing it, and calls
#                               failing only it (those it alone removes)

import sys

import pandas as pd

# Required MAF columns for maftools::read.maf, written first in the MAF outputs
MAF_CORE = ["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
            "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
            "Variant_Type", "Tumor_Sample_Barcode"]

PLACEHOLDER_GENE = "CHIP_NEGATIVE_PLACEHOLDER"


def fail(handle, message):
    """Write the error to the rule log (Snakemake does not capture a script's
    stderr there), then stop."""
    print(f"ERROR: {message}", file=handle, flush=True)
    sys.exit(message)


def write_maf(frame, path):
    columns = MAF_CORE + [c for c in frame.columns if c not in MAF_CORE]
    with open(path, "w") as handle:
        handle.write("#version 2.4\n")
        frame[columns].to_csv(handle, sep="\t", index=False)


def filter_counts(candidates):
    """One row per criterion: calls failing it, and calls failing it and no other."""
    criteria = [c for c in candidates.columns if c.startswith("pass_")]
    rows = []
    for criterion in criteria:
        others = [c for c in criteria if c != criterion]
        failing = ~candidates[criterion]
        only = failing & candidates[others].all(axis=1)
        rows.append({"criterion": criterion, "candidates": len(candidates),
                     "failing": int(failing.sum()), "failing_only_this": int(only.sum())})
    return pd.DataFrame(rows, columns=["criterion", "candidates", "failing", "failing_only_this"])


def main(snakemake, log):
    samples = list(snakemake.params.samples)
    large = snakemake.params.large_clone_vaf

    tables = []
    gnomad_checked = {}
    for sample, path, status_path in zip(samples, snakemake.input.tables, snakemake.input.status):
        table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if table["Tumor_Sample_Barcode"].ne(sample).any():
            fail(log, f"{path}: Tumor_Sample_Barcode differs from sample {sample}")
        tables.append(table)
        status = pd.read_csv(status_path, sep="\t", dtype=str)
        if status.at[0, "sample"] != sample:
            fail(log, f"{status_path}: belongs to {status.at[0, 'sample']}, not {sample}")
        gnomad_checked[sample] = status.at[0, "gnomad_checked"] == "True"
    candidates = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=MAF_CORE)

    flags = [c for c in candidates.columns if c.startswith("pass_")] + \
            ["qualifies", "counts_for_status", "in_pon", "pon_hotspot", "gnomad_checked",
             "is_indel", "in_str", "in_giab_repeat", "in_repeat"]
    for c in flags:
        candidates[c] = candidates[c].astype(str).str.lower().eq("true")
    for c in ["t_vaf", "t_depth", "sor"]:
        candidates[c] = pd.to_numeric(candidates[c])

    candidates.to_csv(snakemake.output.candidates, sep="\t", index=False)

    counts = filter_counts(candidates)
    counts.to_csv(snakemake.output.filter_counts, sep="\t", index=False)
    for row in counts.itertuples(index=False):
        print(f"  {row.criterion}: failing {row.failing}, only this {row.failing_only_this}", file=log)

    counted = candidates[candidates["counts_for_status"]].copy()
    other = candidates[candidates["qualifies"] & ~candidates["counts_for_status"]]
    write_maf(counted, snakemake.output.maf)

    # ---- one row per sample, negatives included ---------------------------------
    rows = []
    for sample in samples:
        calls = counted[counted["Tumor_Sample_Barcode"] == sample]
        others = other[other["Tumor_Sample_Barcode"] == sample]
        rows.append({
            "Tumor_Sample_Barcode": sample,
            "CHIP_status": "CHIP+" if len(calls) else "CHIP-",
            "n_CHIP_mutations": len(calls),
            "max_VAF": round(calls["t_vaf"].max(), 5) if len(calls) else 0.0,
            "sum_VAF": round(calls["t_vaf"].sum(), 5) if len(calls) else 0.0,
            "has_large_clone": bool((calls["t_vaf"] >= large).any()),
            "genes_mutated": ", ".join(sorted(calls["Hugo_Symbol"].unique())),
            "other_genes_qualifying": ", ".join(sorted(others["Hugo_Symbol"].unique())),
            "n_candidates": int((candidates["Tumor_Sample_Barcode"] == sample).sum()),
            "gnomad_checked": gnomad_checked[sample],
        })
    per_sample = pd.DataFrame(rows)
    per_sample.to_csv(snakemake.output.samples, sep="\t", index=False)

    # ---- oncoplot MAF: Silent placeholders keep CHIP- samples in the plot --------
    negatives = per_sample.loc[per_sample["CHIP_status"] == "CHIP-", "Tumor_Sample_Barcode"]
    placeholders = pd.DataFrame({
        "Hugo_Symbol": PLACEHOLDER_GENE, "Chromosome": "chr1", "Start_Position": 1, "End_Position": 1,
        "Reference_Allele": "A", "Tumor_Seq_Allele2": "T", "Variant_Classification": "Silent",
        "Variant_Type": "SNP", "Tumor_Sample_Barcode": list(negatives),
    })
    write_maf(pd.concat([counted, placeholders], ignore_index=True), snakemake.output.oncoplot_maf)

    n_pos = int((per_sample["CHIP_status"] == "CHIP+").sum())
    print(f"samples: {len(per_sample)} | CHIP+: {n_pos} | CHIP-: {len(per_sample) - n_pos}", file=log)
    print(f"candidate calls: {len(candidates)} | qualifying: {int(candidates['qualifies'].sum())} "
          f"| counted for status: {len(counted)}", file=log)
    unchecked = per_sample.loc[~per_sample["gnomad_checked"], "Tumor_Sample_Barcode"].tolist()
    if unchecked:
        print(f"WARNING: gnomAD criterion skipped (no gnomAD column) for: {', '.join(unchecked)}", file=log)


with open(snakemake.log[0], "w") as log_handle:
    main(snakemake, log_handle)
