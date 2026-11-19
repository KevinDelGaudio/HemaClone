# CHIP classification of one sample's ch_final calls (Snakemake script).
#
# Reads the Funcotator MAF and the ch_final VCF it was made from, and writes
# one row per call with every MAF column, the derived values and one
# pass/fail column per criterion. The PoN flag is read from the VCF:
# Funcotator writes one MAF row per alternate allele in VCF order, which is
# checked position by position, so a mismatch stops the step.

import gzip
import math
import re
import sys

import pandas as pd

# GNOMAD_AF: INFO field added by the gnomad_annotate rule (the source used).
# gnomAD_exome_AF / gnomAD_genome_AF: Funcotator's own gnomAD data sources,
# used as well if they were extracted. The largest value counts.
GNOMAD_COLUMN = re.compile(r"^gnomAD(_exome|_genome)?_AF$", re.IGNORECASE)


def log(handle, text):
    print(text, file=handle, flush=True)


def read_maf(path):
    with open(path) as handle:
        lines = [line for line in handle if not line.startswith("#")]
    if not lines:
        sys.exit(f"MAF without a header line: {path}")
    header = lines[0].rstrip("\n").split("\t")
    rows = [line.rstrip("\n").split("\t") for line in lines[1:] if line.strip()]
    for row in rows:
        if len(row) != len(header):
            sys.exit(f"MAF row with {len(row)} fields, header has {len(header)}: {path}")
    return pd.DataFrame(rows, columns=header, dtype=str)


def read_vcf_alleles(path):
    """One entry per alternate allele, in file order: chrom, pos, PON flag.
    Also returns whether the header defines INFO/GNOMAD_AF, i.e. whether the
    gnomAD annotation step ran, independently of how many calls there are."""
    alleles = []
    gnomad_annotated = False
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=GNOMAD_AF,"):
                gnomad_annotated = True
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, alts, info = fields[0], int(fields[1]), fields[4].split(","), fields[7]
            in_pon = "PON" in info.split(";")
            for _ in alts:
                alleles.append((chrom, pos, in_pon))
    return alleles, gnomad_annotated


def binom_p_below_half(alt, depth):
    """One-sided P(X <= alt) for X ~ Binomial(depth, 0.5), as R's binom.test(..., 'less')."""
    if depth <= 0:
        return float("nan")
    return sum(math.comb(depth, i) for i in range(alt + 1)) / 2 ** depth


def gnomad_af(value):
    """Largest number in a Funcotator field ('' or '|'-separated); absent from gnomAD = 0."""
    numbers = []
    for part in re.split(r"[|,]", value or ""):
        try:
            numbers.append(float(part))
        except ValueError:
            pass
    return max(numbers) if numbers else 0.0


def main(snakemake):
    log_handle = open(snakemake.log[0], "w")
    sample = snakemake.wildcards.sample
    chip = snakemake.params.chip
    min_dp = snakemake.params.min_dp
    min_vaf = snakemake.params.min_vaf
    pon_used = snakemake.params.pon_used

    maf = read_maf(snakemake.input.maf)
    alleles, gnomad_annotated = read_vcf_alleles(snakemake.input.vcf)
    log(log_handle, f"{sample}: {len(maf)} MAF rows, {len(alleles)} VCF alleles")

    # ---- MAF rows <-> VCF alleles, checked one by one ---------------------------
    if len(maf) != len(alleles):
        sys.exit(f"{sample}: MAF has {len(maf)} rows but the VCF has {len(alleles)} alleles")
    for i, (chrom, pos, _) in enumerate(alleles):
        maf_chrom, maf_start = maf.at[i, "Chromosome"], int(maf.at[i, "Start_Position"])
        if maf_chrom != chrom or maf_start not in (pos, pos + 1):
            sys.exit(f"{sample}: MAF row {i + 1} ({maf_chrom}:{maf_start}) does not match "
                     f"VCF allele {i + 1} ({chrom}:{pos})")

    required = ["Hugo_Symbol", "Variant_Classification", "Protein_Change", "t_alt_count", "t_ref_count"]
    missing = [c for c in required if c not in maf.columns]
    if missing:
        sys.exit(f"{sample}: MAF lacks columns {missing}")

    # gnomAD availability: the annotated VCF header says whether the annotation
    # ran. Funcotator writes no INFO columns for a sample without calls, so an
    # empty MAF lacking GNOMAD_AF is expected; a MAF with calls lacking it is not.
    gnomad_cols = [c for c in maf.columns if GNOMAD_COLUMN.match(c)]
    if gnomad_annotated and len(maf) > 0 and not gnomad_cols:
        sys.exit(f"{sample}: the VCF carries INFO/GNOMAD_AF but the MAF has no gnomAD column")
    gnomad_checked = gnomad_annotated or bool(gnomad_cols)
    if gnomad_cols:
        log(log_handle, f"gnomAD columns: {', '.join(gnomad_cols)}")
    elif gnomad_checked:
        log(log_handle, "gnomAD annotation present; no calls to check")
    else:
        log(log_handle, "WARNING: no gnomAD annotation (GNOMAD_AF missing from the VCF header: was the "
                        "gnomad_annotate step skipped?); the gnomAD criterion is skipped for this sample")

    # ---- derived values -----------------------------------------------------------
    maf["Tumor_Sample_Barcode"] = sample
    alt = pd.to_numeric(maf["t_alt_count"], errors="coerce").fillna(0).astype(int)
    ref = pd.to_numeric(maf["t_ref_count"], errors="coerce").fillna(0).astype(int)
    depth = alt + ref
    vaf = (alt / depth.where(depth > 0)).fillna(0.0)

    maf["t_depth"] = depth
    maf["t_vaf"] = vaf.round(5)
    maf["binom_p_below_half"] = [binom_p_below_half(a, d) for a, d in zip(alt, depth)]
    if gnomad_cols:
        maf["gnomad_af"] = [max(gnomad_af(v) for v in values)
                            for values in maf[gnomad_cols].itertuples(index=False)]
    else:
        maf["gnomad_af"] = float("nan")
    maf["in_pon"] = [in_pon for _, _, in_pon in alleles] if pon_used else False

    hotspots = [h.split(":", 1) for h in chip["pon_exempt_hotspots"]]
    maf["pon_hotspot"] = [
        any(gene == g and str(change).startswith(prefix) for g, prefix in hotspots)
        for gene, change in zip(maf["Hugo_Symbol"], maf["Protein_Change"])
    ]

    # ---- criteria: one boolean column each -------------------------------------
    maf["pass_class"] = maf["Variant_Classification"].isin(chip["variant_classes"])
    maf["pass_alt_reads"] = alt >= chip["min_alt_reads"]
    maf["pass_depth"] = depth >= min_dp
    maf["pass_min_vaf"] = vaf >= min_vaf
    maf["pass_max_vaf"] = vaf <= chip["max_vaf"]
    maf["pass_not_germline"] = maf["binom_p_below_half"] < chip["max_binom_p"]
    maf["pass_gnomad"] = (maf["gnomad_af"] <= chip["max_gnomad_af"]) if gnomad_cols else True
    maf["pass_pon"] = ~maf["in_pon"] | maf["pon_hotspot"]

    criteria = [c for c in maf.columns if c.startswith("pass_")]
    maf["qualifies"] = maf[criteria].all(axis=1)
    maf["counts_for_status"] = maf["qualifies"] & ~maf["Hugo_Symbol"].isin(chip["status_exclude_genes"])
    maf["gnomad_checked"] = gnomad_checked

    maf.to_csv(snakemake.output.candidates, sep="\t", index=False)

    # Per-sample status, independent of how many calls the sample has
    with open(snakemake.output.status, "w") as handle:
        handle.write("sample\tgnomad_checked\tgnomad_columns\tn_candidates\n")
        handle.write(f"{sample}\t{gnomad_checked}\t{','.join(gnomad_cols) or '-'}\t{len(maf)}\n")

    log(log_handle, f"qualifying calls: {int(maf['qualifies'].sum())} "
                    f"(counted for CHIP status: {int(maf['counts_for_status'].sum())})")
    for c in criteria:
        log(log_handle, f"  failing {c}: {int((~maf[c]).sum())}")
    log_handle.close()


main(snakemake)
