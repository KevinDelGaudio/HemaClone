#!/usr/bin/env python3
"""
Download and build the hg38 resources HemaClone needs, into one folder whose
file names match the defaults in config/config.yaml.

Normally launched through Docker by hemaclone.sh:
  bash hemaclone.sh resources --resources-dir /data/hg38_resources \
      --targets /data/kits/exome_targets.hg38.bed \
      --cores 16 --mem-gb 96 --tmpdir /scratch/tmp

Steps, each skipped when its output is already recorded in MANIFEST.tsv
with the same size and checksum:
  1. GATK reference bundle: FASTA, index, dictionary, ALT contig list, dbSNP,
     Mills indels, WGS calling regions
  2. GATK somatic resources: gnomAD allele frequencies, common sites, panel
     of normals
  3. GIAB tandem repeats and homopolymers BED, bgzipped and tabix-indexed
  4. somalier sites (sex and sample identity check)
  5. bwa-mem2 index of the FASTA
  6. Funcotator somatic data sources
  7. Exome intervals: the kit BED given with --targets, as an interval list
  8. CH gene panel BED: gene spans from the Funcotator GENCODE annotation,
     intersected with the kit targets and padded

Downloads resume after an interruption and are checked against the size and,
when the server provides one, the MD5 checksum. MANIFEST.tsv records every
file with its source, size, MD5 and date.

Settings come from config.yaml, --settings FILE and --set KEY=VALUE, as for a
run; the resources section holds the keys used only here.
"""

import argparse
import base64
import datetime
import glob
import hashlib
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_batches as rb  # noqa: E402  (settings handling shared with runs)

MANIFEST_NAME = "MANIFEST.tsv"
MANIFEST_FIELDS = ["file", "source", "bytes", "md5", "date"]
CHUNK = 8 * 1024 * 1024

# (file name in the bundle, file name saved). The ALT contig list is saved
# under the name bwa-mem2 looks for: <FASTA>.alt, next to the index.
BUNDLE_FILES = [
    ("Homo_sapiens_assembly38.fasta", "Homo_sapiens_assembly38.fasta"),
    ("Homo_sapiens_assembly38.fasta.fai", "Homo_sapiens_assembly38.fasta.fai"),
    ("Homo_sapiens_assembly38.dict", "Homo_sapiens_assembly38.dict"),
    ("Homo_sapiens_assembly38.fasta.64.alt", "Homo_sapiens_assembly38.fasta.alt"),
    ("Homo_sapiens_assembly38.dbsnp138.vcf", "Homo_sapiens_assembly38.dbsnp138.vcf"),
    ("Homo_sapiens_assembly38.dbsnp138.vcf.idx", "Homo_sapiens_assembly38.dbsnp138.vcf.idx"),
    ("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"),
    ("Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"),
    ("wgs_calling_regions.hg38.interval_list", "wgs_calling_regions.hg38.interval_list"),
]

SOMATIC_FILES = [
    ("af-only-gnomad.hg38.vcf.gz", "af-only-gnomad.hg38.vcf.gz"),
    ("af-only-gnomad.hg38.vcf.gz.tbi", "af-only-gnomad.hg38.vcf.gz.tbi"),
    ("small_exac_common_3.hg38.vcf.gz", "small_exac_common_3.hg38.vcf.gz"),
    ("small_exac_common_3.hg38.vcf.gz.tbi", "small_exac_common_3.hg38.vcf.gz.tbi"),
    ("1000g_pon.hg38.vcf.gz", "1000g_pon.hg38.vcf.gz"),
    ("1000g_pon.hg38.vcf.gz.tbi", "1000g_pon.hg38.vcf.gz.tbi"),
]

FASTA = "Homo_sapiens_assembly38.fasta"
DICT = "Homo_sapiens_assembly38.dict"
BWA_INDEX_FILE = FASTA + ".bwt.2bit.64"

# Keys with a dedicated flag here
KEYS_WITH_FLAGS = {
    "reference.dir": "--resources-dir",
    "resources.targets": "--targets",
    "tmpdir": "--tmpdir",
    "machine.cores": "--cores",
    "machine.mem_gb": "--mem-gb",
}


def fail(message):
    sys.exit(f"ERROR: {message}")


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def step(text):
    print(f"\n=== {text}", flush=True)


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Download and build the hg38 resources HemaClone needs.")
    parser.add_argument("--resources-dir", required=True, help="folder to fill (created if missing)")
    parser.add_argument("--targets", default=None, help="exome kit target BED, hg38 (resources.targets)")
    parser.add_argument("--cores", type=int, required=True, help="CPU limit of the container")
    parser.add_argument("--mem-gb", type=float, required=True, help="memory limit, GB")
    parser.add_argument("--tmpdir", default=None, help="scratch folder (default: tmpdir)")
    parser.add_argument("--settings", default=None, help="YAML file with config.yaml keys")
    parser.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                        help="override one config.yaml key; repeatable")
    return parser.parse_args(argv)


def load_config(args):
    with open(rb.CONFIG) as handle:
        config = yaml.safe_load(handle)
    if args.settings:
        rb.deep_update(config, rb.load_settings_file(args.settings, config))
    for item in args.set:
        key = item.split("=", 1)[0]
        if key in KEYS_WITH_FLAGS:
            fail(f"use {KEYS_WITH_FLAGS[key]} instead of --set {key}")
    rb.deep_update(config, rb.parse_set_overrides(args.set, config))
    if args.targets:
        config["resources"]["targets"] = args.targets
    if args.tmpdir:
        config["tmpdir"] = args.tmpdir
    return config


# -----------------------------------------------------------------------------
# Manifest: one row per file, rewritten after every step
# -----------------------------------------------------------------------------
def read_manifest(dest):
    path = os.path.join(dest, MANIFEST_NAME)
    rows = {}
    if os.path.isfile(path):
        with open(path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            if header != MANIFEST_FIELDS:
                fail(f"{path}: unexpected header {header}, expected {MANIFEST_FIELDS}")
            for line in handle:
                values = line.rstrip("\n").split("\t")
                rows[values[0]] = dict(zip(MANIFEST_FIELDS, values))
    return rows


def write_manifest(dest, rows):
    path = os.path.join(dest, MANIFEST_NAME)
    with open(f"{path}.part", "w") as handle:
        handle.write("\t".join(MANIFEST_FIELDS) + "\n")
        for name in sorted(rows):
            handle.write("\t".join(str(rows[name][f]) for f in MANIFEST_FIELDS) + "\n")
    os.replace(f"{path}.part", path)


def md5_of(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def recorded(rows, dest, name):
    """True when the file exists with the size recorded in the manifest."""
    path = os.path.join(dest, name)
    row = rows.get(name)
    return row is not None and os.path.isfile(path) and os.path.getsize(path) == int(row["bytes"])


def record_file(rows, dest, name, source):
    path = os.path.join(dest, name)
    rows[name] = {"file": name, "source": source, "bytes": os.path.getsize(path),
                  "md5": md5_of(path), "date": now()}


# -----------------------------------------------------------------------------
# Downloads
# -----------------------------------------------------------------------------
def remote_info(url):
    """Size and MD5 (hex, or None) from a HEAD request. Google Cloud Storage
    gives the MD5 in x-goog-hash for files uploaded in one piece."""
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            size = int(response.headers["Content-Length"])
            md5 = None
            for part in response.headers.get_all("x-goog-hash") or []:
                for item in part.split(","):
                    item = item.strip()
                    if item.startswith("md5="):
                        md5 = base64.b64decode(item[4:]).hex()
            return size, md5
    except (urllib.error.URLError, KeyError, ValueError) as error:
        fail(f"cannot reach {url}: {error}")


def download(url, path, size):
    """Download to PATH.part, resuming a partial file when the server allows it."""
    part = f"{path}.part"
    done = os.path.getsize(part) if os.path.isfile(part) else 0
    if done > size:
        os.remove(part)
        done = 0
    # Interrupted after the last byte, before the rename: nothing left to fetch
    # (a range request past the end would be refused)
    if done == size:
        os.replace(part, path)
        return
    request = urllib.request.Request(url)
    if done:
        request.add_header("Range", f"bytes={done}-")
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            if done and response.status != 206:
                # The server sent the whole file: start again
                done = 0
            with open(part, "ab" if done else "wb") as out:
                copied = done
                last_report = 0
                for block in iter(lambda: response.read(CHUNK), b""):
                    out.write(block)
                    copied += len(block)
                    if copied - last_report >= 1024 ** 3:
                        print(f"    {copied / 1024 ** 3:.1f} of {size / 1024 ** 3:.1f} GB", flush=True)
                        last_report = copied
    except (urllib.error.URLError, OSError) as error:
        fail(f"download of {url} interrupted: {error}. Re-run the same command: it resumes "
             f"from {part}")
    os.replace(part, path)


def fetch_files(base_url, files, dest, rows):
    for remote_name, local_name in files:
        url = f"{base_url.rstrip('/')}/{remote_name}"
        size, md5 = remote_info(url)
        path = os.path.join(dest, local_name)
        if recorded(rows, dest, local_name) and int(rows[local_name]["bytes"]) == size \
                and (md5 is None or rows[local_name]["md5"] == md5):
            print(f"  present   {local_name}")
            continue
        print(f"  download  {local_name} ({size / 1024 ** 2:,.0f} MB)", flush=True)
        download(url, path, size)
        if os.path.getsize(path) != size:
            fail(f"{local_name}: downloaded {os.path.getsize(path)} bytes, expected {size}")
        record_file(rows, dest, local_name, url)
        if md5 is not None and rows[local_name]["md5"] != md5:
            got = rows.pop(local_name)["md5"]
            os.remove(path)
            write_manifest(dest, rows)
            fail(f"{local_name}: MD5 {got} does not match the server's {md5}; "
                 "the file was removed, re-run to download it again")
        write_manifest(dest, rows)


# -----------------------------------------------------------------------------
# Tools
# -----------------------------------------------------------------------------
def run_tool(command, log_path):
    """Run a tool with its output in LOG_PATH; stop with the log's end on failure."""
    print(f"  $ {' '.join(command)}", flush=True)
    with open(log_path, "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        with open(log_path) as log:
            tail = log.readlines()[-20:]
        hint = ""
        if result.returncode in (137, -9):
            hint = " It was killed, most likely for lack of memory: raise --mem-gb."
        fail(f"{command[0]} failed (exit {result.returncode}).{hint} Last lines of {log_path}:\n"
             + "".join(tail))


def get_repeats_bed(url, dest, rows, name, logs):
    """GIAB stratification BED: the download is kept under giab/, the copy
    the workflow reads is bgzipped and indexed with tabix."""
    source = os.path.join("giab", os.path.basename(url))
    os.makedirs(os.path.join(dest, "giab"), exist_ok=True)
    size, md5 = remote_info(url)
    if recorded(rows, dest, source) and int(rows[source]["bytes"]) == size \
            and recorded(rows, dest, name) and recorded(rows, dest, f"{name}.tbi"):
        print(f"  present   {name}")
        return
    print(f"  download  {source} ({size / 1024 ** 2:,.0f} MB)", flush=True)
    download(url, os.path.join(dest, source), size)
    record_file(rows, dest, source, url)
    if md5 is not None and rows[source]["md5"] != md5:
        got = rows.pop(source)["md5"]
        os.remove(os.path.join(dest, source))
        write_manifest(dest, rows)
        fail(f"{source}: MD5 {got} does not match the server's {md5}; re-run to download it again")
    target = os.path.join(dest, name)
    log_path = os.path.join(logs, "repeats_bed.log")
    print(f"  $ gzip -dc {source} | bgzip -c > {name}; tabix -p bed {name}", flush=True)
    with open(log_path, "w") as log, open(f"{target}.part", "wb") as out:
        unzip = subprocess.Popen(["gzip", "-dc", os.path.join(dest, source)],
                                 stdout=subprocess.PIPE, stderr=log)
        zipped = subprocess.run(["bgzip", "-c"], stdin=unzip.stdout, stdout=out, stderr=log)
        unzip.stdout.close()
        if unzip.wait() != 0 or zipped.returncode != 0:
            out.close()
            os.remove(f"{target}.part")
            fail(f"recompressing {source} failed; see {log_path}")
    os.replace(f"{target}.part", target)
    run_tool(["tabix", "-f", "-p", "bed", target], log_path)
    record_file(rows, dest, name, f"built: bgzip of {source}")
    record_file(rows, dest, f"{name}.tbi", f"built: tabix -p bed {name}")
    write_manifest(dest, rows)


def build_bwa_index(dest, rows, logs):
    fasta = os.path.join(dest, FASTA)
    index_files = [f"{FASTA}{suffix}" for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]
    if all(recorded(rows, dest, name) for name in index_files):
        print("  present   bwa-mem2 index")
        return
    run_tool(["bwa-mem2", "index", fasta], os.path.join(logs, "bwa_mem2_index.log"))
    for name in index_files:
        record_file(rows, dest, name, f"built: bwa-mem2 index {FASTA}")
    write_manifest(dest, rows)


def get_funcotator(dest, rows, folder_name, tmpdir, logs):
    marker = f"{folder_name}/VERSION"
    if recorded(rows, dest, marker) and os.path.isdir(os.path.join(dest, folder_name)):
        print(f"  present   {folder_name} ({rows[marker]['source']})")
        return
    work = os.path.join(dest, ".funcotator_download")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    tarball = os.path.join(work, "funcotator_dataSources.tar.gz")
    run_tool(["gatk", "--java-options", f"-Djava.io.tmpdir={tmpdir}",
              "FuncotatorDataSourceDownloader", "--somatic", "--hg38",
              "--validate-integrity", "--extract-after-download", "--output", tarball],
             os.path.join(logs, "funcotator_download.log"))
    extracted = [d for d in glob.glob(os.path.join(work, "funcotator_dataSources*")) if os.path.isdir(d)]
    if len(extracted) != 1:
        fail(f"expected one extracted Funcotator folder in {work}, found {extracted}")
    release = os.path.basename(extracted[0])
    target = os.path.join(dest, folder_name)
    shutil.rmtree(target, ignore_errors=True)
    os.replace(extracted[0], target)
    shutil.rmtree(work)
    with open(os.path.join(target, "VERSION"), "w") as handle:
        handle.write(release + "\n")
    record_file(rows, dest, marker, f"gatk FuncotatorDataSourceDownloader --somatic --hg38: {release}")
    write_manifest(dest, rows)


# -----------------------------------------------------------------------------
# Intervals
# -----------------------------------------------------------------------------
def contig_lengths(dest):
    """Contig name -> length, in dictionary order (from the FASTA index)."""
    lengths = {}
    with open(os.path.join(dest, f"{FASTA}.fai")) as handle:
        for line in handle:
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)
    return lengths


def read_bed(path, lengths):
    """BED intervals (0-based, half-open) as (contig, start, end); stops on
    contigs missing from the reference or coordinates past a contig's end."""
    intervals = []
    unknown = set()
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                fail(f"{path} line {number}: expected at least 3 tab-separated columns")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            if contig not in lengths:
                unknown.add(contig)
                continue
            if start < 0 or end > lengths[contig] or start >= end:
                fail(f"{path} line {number}: {contig}:{start}-{end} is not a valid interval "
                     f"(contig length {lengths[contig]})")
            intervals.append((contig, start, end))
    if unknown:
        fail(f"{path}: contigs not in the reference: {', '.join(sorted(unknown)[:10])}. "
             "The BED must use hg38 names with the chr prefix (chr1, chrX).")
    if not intervals:
        fail(f"{path}: no intervals")
    print(f"  {path}: {len(intervals)} intervals")
    return intervals


def merge(intervals):
    """Merge overlapping or touching intervals on the same contig (input sorted)."""
    merged = []
    for contig, start, end in intervals:
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1] = (contig, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((contig, start, end))
    return merged


def gene_spans(gtf, genes, lengths):
    """Gene name -> (contig, start, end) from the 'gene' records of a GENCODE
    GTF. Stops when a gene is missing or found at more than one locus."""
    wanted = set(genes)
    found = {}
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            name = None
            for attribute in fields[8].split(";"):
                attribute = attribute.strip()
                if attribute.startswith("gene_name "):
                    name = attribute.split(" ", 1)[1].strip('"')
                    break
            if name in wanted and fields[0] in lengths:
                found.setdefault(name, set()).add((fields[0], int(fields[3]) - 1, int(fields[4])))
    missing = sorted(wanted - set(found))
    if missing:
        fail(f"genes not found in {gtf}: {', '.join(missing)}")
    ambiguous = {name: loci for name, loci in found.items() if len(loci) > 1}
    if ambiguous:
        listing = "; ".join(f"{name}: {sorted(loci)}" for name, loci in sorted(ambiguous.items()))
        fail(f"genes at more than one locus in {gtf}: {listing}")
    return {name: next(iter(loci)) for name, loci in found.items()}


def build_panel(dest, rows, config, targets, logs):
    resources = config["resources"]
    panel_name = config["reference"]["gene_panel_bed"]
    genes = resources["panel_genes"]
    padding = resources["panel_padding"]
    funcotator = os.path.join(dest, config["reference"]["funcotator_ds"])
    gtfs = sorted(glob.glob(os.path.join(funcotator, "gencode", "hg38", "*.gtf")))
    if len(gtfs) != 1:
        fail(f"expected one GENCODE GTF in {funcotator}/gencode/hg38, found {gtfs}")
    gtf = gtfs[0]
    lengths = contig_lengths(dest)
    order = {name: i for i, name in enumerate(lengths)}

    spans = gene_spans(gtf, genes, lengths)
    target_intervals = sorted(read_bed(targets, lengths), key=lambda t: (order[t[0]], t[1]))

    panel = []          # (contig, start, end, gene)
    uncovered = []
    for gene in genes:
        contig, g_start, g_end = spans[gene]
        pieces = []
        for t_contig, t_start, t_end in target_intervals:
            if t_contig == contig and t_start < g_end and t_end > g_start:
                start = max(0, max(t_start, g_start) - padding)
                end = min(lengths[contig], min(t_end, g_end) + padding)
                pieces.append((contig, start, end))
        if not pieces:
            uncovered.append(gene)
            continue
        for contig_, start, end in merge(sorted(pieces, key=lambda p: p[1])):
            panel.append((contig_, start, end, gene))
    panel.sort(key=lambda p: (order[p[0]], p[1], p[3]))

    path = os.path.join(dest, panel_name)
    with open(f"{path}.part", "w") as handle:
        for contig, start, end, gene in panel:
            handle.write(f"{contig}\t{start}\t{end}\t{gene}\n")
    os.replace(f"{path}.part", path)

    covered = len(genes) - len(uncovered)
    total_bp = sum(end - start for _, start, end, _ in panel)
    print(f"  {panel_name}: {len(panel)} intervals, {total_bp:,} bp, "
          f"{covered} of {len(genes)} genes covered by the kit (padding {padding} bp)")
    if uncovered:
        print(f"  WARNING: genes with no kit target, not called: {', '.join(uncovered)}", flush=True)
    record_file(rows, dest, panel_name,
                f"built: {len(genes)} genes from {os.path.basename(gtf)} x {os.path.basename(targets)}, "
                f"padding {padding}; not covered: {', '.join(uncovered) or 'none'}")
    write_manifest(dest, rows)


def build_exome_intervals(dest, rows, config, targets, logs):
    name = config["reference"]["exome_intervals"]
    # Check names and coordinates first, for a clearer error than the tool's
    read_bed(targets, contig_lengths(dest))
    kit_copy = os.path.join("targets", os.path.basename(targets))
    os.makedirs(os.path.join(dest, "targets"), exist_ok=True)
    if os.path.abspath(targets) != os.path.abspath(os.path.join(dest, kit_copy)):
        shutil.copyfile(targets, os.path.join(dest, kit_copy))
    record_file(rows, dest, kit_copy, f"copied from {targets}")
    run_tool(["gatk", "BedToIntervalList",
              "--INPUT", os.path.join(dest, kit_copy),
              "--SEQUENCE_DICTIONARY", os.path.join(dest, DICT),
              "--OUTPUT", os.path.join(dest, name)],
             os.path.join(logs, "bed_to_interval_list.log"))
    record_file(rows, dest, name, f"built: gatk BedToIntervalList {kit_copy}")
    write_manifest(dest, rows)
    return os.path.join(dest, kit_copy)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main(argv):
    args = parse_args(argv)
    config = load_config(args)
    resources = config["resources"]
    dest = os.path.abspath(args.resources_dir)
    tmpdir = os.path.abspath(config["tmpdir"])
    targets = resources["targets"]

    problems = []
    if not targets:
        problems.append("no exome kit targets: pass --targets KIT.bed (hg38, chr names)")
    elif not os.path.isfile(targets):
        problems.append(f"targets BED not found: {targets}")
    for tool in ("gatk", "bwa-mem2", "bgzip", "tabix"):
        if shutil.which(tool) is None:
            problems.append(f"{tool} not found on PATH")
    if not isinstance(resources["panel_padding"], int) or resources["panel_padding"] < 0:
        problems.append(f"resources.panel_padding must be a whole number >= 0, got {resources['panel_padding']!r}")
    if not resources["panel_genes"]:
        problems.append("resources.panel_genes is empty")
    if problems:
        fail("\n  " + "\n  ".join(problems))

    os.makedirs(dest, exist_ok=True)
    os.makedirs(tmpdir, exist_ok=True)
    logs = os.path.join(dest, "logs")
    os.makedirs(logs, exist_ok=True)
    if not os.access(dest, os.W_OK):
        fail(f"resources folder not writable: {dest}")

    print("Resources")
    print(f"  folder       : {dest}")
    print(f"  kit targets  : {targets}")
    print(f"  panel        : {len(resources['panel_genes'])} genes, padding {resources['panel_padding']} bp")
    print(f"  bundle       : {resources['bundle_url']}")
    print(f"  somatic      : {resources['somatic_url']}")
    print(f"  repeats      : {resources['repeats_url']}")
    print(f"  somalier     : {resources['somalier_sites_url']}")
    print(f"  tmpdir       : {tmpdir}")

    rows = read_manifest(dest)

    step("1/8 GATK reference bundle")
    fetch_files(resources["bundle_url"], BUNDLE_FILES, dest, rows)
    step("2/8 GATK somatic resources")
    fetch_files(resources["somatic_url"], SOMATIC_FILES, dest, rows)
    step("3/8 GIAB tandem repeats and homopolymers")
    get_repeats_bed(resources["repeats_url"], dest, rows, config["reference"]["repeats_bed"], logs)
    step("4/8 somalier sites")
    sites_url = resources["somalier_sites_url"]
    fetch_files(os.path.dirname(sites_url),
                [(os.path.basename(sites_url), config["reference"]["somalier_sites"])], dest, rows)
    step("5/8 bwa-mem2 index (slow, and needs much memory for hg38)")
    build_bwa_index(dest, rows, logs)
    step("6/8 Funcotator somatic data sources")
    get_funcotator(dest, rows, config["reference"]["funcotator_ds"], tmpdir, logs)
    step("7/8 Exome intervals")
    kit = build_exome_intervals(dest, rows, config, targets, logs)
    step("8/8 CH gene panel")
    build_panel(dest, rows, config, kit, logs)

    defaults = {"fasta": FASTA, "dbsnp": "Homo_sapiens_assembly38.dbsnp138.vcf",
                "mills": "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "gnomad": "af-only-gnomad.hg38.vcf.gz",
                "common_variants": "small_exac_common_3.hg38.vcf.gz"}
    differing = [f"--set reference.{key}={value}" for key, value in defaults.items()
                 if config["reference"][key] != value]

    print(f"\nDone: {len(rows)} files recorded in {os.path.join(dest, MANIFEST_NAME)}")
    print(f"Run with: --resources-dir {dest}")
    if differing:
        print("Your settings use other file names than the ones downloaded; add: " + " ".join(differing))


if __name__ == "__main__":
    main(sys.argv[1:])
