# CHIP classification of one sample's ch_final calls (Snakemake script).
#
# Reads the Funcotator MAF and the ch_final VCF it was made from, and writes
# one row per call with every MAF column, the derived values and one
# pass/fail column per criterion. The PoN flag, the strand counts (FORMAT/SB)
# and the repeat flag (INFO/STR) are read from the VCF: Funcotator writes one
# MAF row per alternate allele in VCF order, which is checked position by
# position, so a mismatch stops the step. Indels are also looked up in the
# GIAB tandem-repeat and homopolymer BED.

import gzip
import math
import re
import subprocess
import sys

import pandas as pd

# GNOMAD_AF: INFO field added by the gnomad_annotate rule (the source used).
# gnomAD_exome_AF / gnomAD_genome_AF: Funcotator's own gnomAD data sources,
# used as well if they were extracted. The largest value counts.
GNOMAD_COLUMN = re.compile(r"^gnomAD(_exome|_genome)?_AF$", re.IGNORECASE)


def log(handle, text):
    print(text, file=handle, flush=True)


def fail(handle, message):
    """Write the error to the rule log (Snakemake does not capture a script's
    stderr there), then stop."""
    log(handle, f"ERROR: {message}")
    sys.exit(message)


def read_maf(path, log_handle):
    with open(path) as handle:
        lines = [line for line in handle if not line.startswith("#")]
    if not lines:
        fail(log_handle, f"MAF without a header line: {path}")
    header = lines[0].rstrip("\n").split("\t")
    rows = [line.rstrip("\n").split("\t") for line in lines[1:] if line.strip()]
    for row in rows:
        if len(row) != len(header):
            fail(log_handle, f"MAF row with {len(row)} fields, header has {len(header)}: {path}")
    return pd.DataFrame(rows, columns=header, dtype=str)


def read_vcf_alleles(path):
    """One entry per alternate allele, in file order, with chrom, pos, ref,
    alt, the PON and STR flags and the sample's FORMAT/SB counts (one set per
    site, shared by its alleles). Also returns which header lines are present:
    INFO/GNOMAD_AF (the gnomAD annotation step ran), FORMAT/SB and INFO/STR."""
    alleles = []
    header = {"gnomad": False, "sb": False, "str": False}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=GNOMAD_AF,"):
                header["gnomad"] = True
            elif line.startswith("##FORMAT=<ID=SB,"):
                header["sb"] = True
            elif line.startswith("##INFO=<ID=STR,"):
                header["str"] = True
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, info = fields[0], int(fields[1]), fields[3], fields[4].split(","), fields[7]
            flags = info.split(";")
            sb = None
            if len(fields) > 9:
                values = dict(zip(fields[8].split(":"), fields[9].split(":")))
                if values.get("SB", ".") != ".":
                    sb = [int(x) for x in values["SB"].split(",")]
            for alt in alts:
                alleles.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                                "in_pon": "PON" in flags, "in_str": "STR" in flags, "sb": sb})
    return alleles, header


def binom_p_below_half(alt, depth):
    """One-sided P(X <= alt) for X ~ Binomial(depth, 0.5), as R's binom.test(..., 'less')."""
    if depth <= 0:
        return float("nan")
    return sum(math.comb(depth, i) for i in range(alt + 1)) / 2 ** depth


def strand_odds_ratio(ref_fwd, ref_rev, alt_fwd, alt_rev):
    """GATK's StrandOddsRatio, with 1 added to each count."""
    a, b, c, d = ref_fwd + 1, ref_rev + 1, alt_fwd + 1, alt_rev + 1
    ratio = (a * d) / (b * c)
    symmetric = ratio + 1 / ratio
    ref_ratio = min(a, b) / max(a, b)
    alt_ratio = min(c, d) / max(c, d)
    return math.log(symmetric) + math.log(ref_ratio) - math.log(alt_ratio)


def repeat_overlaps(bed, regions):
    """For each (chrom, start, end) region, 1-based inclusive, whether it
    overlaps an interval of the tabix-indexed BED."""
    if not regions:
        return []
    queries = [f"{chrom}:{start}-{end}" for chrom, start, end in regions]
    result = subprocess.run(["tabix", bed, *queries], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"tabix {bed} failed: {result.stderr.strip()}")
    intervals = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        intervals.append((fields[0], int(fields[1]), int(fields[2])))
    # BED is 0-based half-open: [s, e) covers 1-based positions s+1 .. e
    return [any(c == chrom and s < end and e >= start for c, s, e in intervals)
            for chrom, start, end in regions]


def gnomad_af(value):
    """Largest number in a Funcotator field ('' or '|'-separated); absent from gnomAD = 0."""
    numbers = []
    for part in re.split(r"[|,]", value or ""):
        try:
            numbers.append(float(part))
        except ValueError:
            pass
    return max(numbers) if numbers else 0.0


def main(snakemake, log_handle):
    sample = snakemake.wildcards.sample
    chip = snakemake.params.chip
    min_dp = snakemake.params.min_dp
    min_vaf = snakemake.params.min_vaf
    pon_used = snakemake.params.pon_used

    maf = read_maf(snakemake.input.maf, log_handle)
    alleles, header = read_vcf_alleles(snakemake.input.vcf)
    gnomad_annotated = header["gnomad"]
    log(log_handle, f"{sample}: {len(maf)} MAF rows, {len(alleles)} VCF alleles")

    # ---- MAF rows <-> VCF alleles, checked one by one ---------------------------
    if len(maf) != len(alleles):
        fail(log_handle, f"{sample}: MAF has {len(maf)} rows but the VCF has {len(alleles)} alleles")
    for i, allele in enumerate(alleles):
        chrom, pos = allele["chrom"], allele["pos"]
        maf_chrom, maf_start = maf.at[i, "Chromosome"], int(maf.at[i, "Start_Position"])
        if maf_chrom != chrom or maf_start not in (pos, pos + 1):
            fail(log_handle, f"{sample}: MAF row {i + 1} ({maf_chrom}:{maf_start}) does not match "
                     f"VCF allele {i + 1} ({chrom}:{pos})")

    required = ["Hugo_Symbol", "Variant_Classification", "Protein_Change", "t_alt_count", "t_ref_count"]
    missing = [c for c in required if c not in maf.columns]
    if missing:
        fail(log_handle, f"{sample}: MAF lacks columns {missing}")

    # gnomAD availability: the annotated VCF header says whether the annotation
    # ran. Funcotator writes no INFO columns for a sample without calls, so an
    # empty MAF lacking GNOMAD_AF is expected; a MAF with calls lacking it is not.
    gnomad_cols = [c for c in maf.columns if GNOMAD_COLUMN.match(c)]
    if gnomad_annotated and len(maf) > 0 and not gnomad_cols:
        fail(log_handle, f"{sample}: the VCF carries INFO/GNOMAD_AF but the MAF has no gnomAD column")
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
    maf["in_pon"] = [a["in_pon"] for a in alleles] if pon_used else False

    # Strand counts: Mutect2 writes FORMAT/SB for every call
    if alleles and not header["sb"]:
        fail(log_handle, f"{sample}: the VCF has no FORMAT/SB (Mutect2 strand counts); "
                         "the strand-bias criterion cannot be applied")
    missing_sb = [f"{a['chrom']}:{a['pos']}" for a in alleles if a["sb"] is None or len(a["sb"]) != 4]
    if missing_sb:
        fail(log_handle, f"{sample}: FORMAT/SB missing or not 4 counts at {', '.join(missing_sb[:5])}")
    for i, name in enumerate(["sb_ref_fwd", "sb_ref_rev", "sb_alt_fwd", "sb_alt_rev"]):
        maf[name] = [a["sb"][i] for a in alleles]
    maf["sor"] = [round(strand_odds_ratio(*a["sb"]), 4) for a in alleles]

    # Indels in tandem repeats or homopolymers: Mutect2 INFO/STR or the GIAB BED
    if alleles and not header["str"]:
        log(log_handle, "INFO/STR absent from the VCF header: repeats taken from the GIAB BED only")
    maf["is_indel"] = [len(a["ref"]) != len(a["alt"]) for a in alleles]
    maf["in_str"] = [a["in_str"] for a in alleles]
    indel_rows = [i for i, a in enumerate(alleles) if len(a["ref"]) != len(a["alt"])]
    # The indel's span on the reference: the anchor base, the bases it deletes
    # and the next base, where an insertion goes
    regions = [(alleles[i]["chrom"], alleles[i]["pos"], alleles[i]["pos"] + len(alleles[i]["ref"]))
               for i in indel_rows]
    try:
        overlaps = repeat_overlaps(snakemake.input.repeats, regions)
    except RuntimeError as error:
        fail(log_handle, f"{sample}: {error}")
    in_giab = [False] * len(alleles)
    for i, overlap in zip(indel_rows, overlaps):
        in_giab[i] = overlap
    maf["in_giab_repeat"] = in_giab
    maf["in_repeat"] = maf["is_indel"] & (maf["in_str"] | maf["in_giab_repeat"])

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
    maf["pass_strand"] = maf["sor"] <= chip["max_sor"]
    maf["pass_indel"] = ~maf["is_indel"] | ((vaf >= chip["indel_min_vaf"])
                                            & (alt >= chip["indel_min_alt_reads"]))
    maf["pass_repeat_indel"] = ~maf["in_repeat"] | ((vaf >= chip["repeat_indel_min_vaf"])
                                                    & (alt >= chip["repeat_indel_min_alt_reads"]))

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


with open(snakemake.log[0], "w") as log_handle:
    main(snakemake, log_handle)
