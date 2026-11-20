#!/usr/bin/env python3
"""
Run HemaClone in batches of N samples, with every setting available
as a launch flag.

Each batch runs to completion (alignment through MAFs and per-sample QC)
before the next batch starts, so no more than N samples are ever in progress.
After the last batch, the cohort QC report (MultiQC, qc_flags.tsv) is built
once over the whole sample sheet.

Finished samples are recognised by Snakemake and skipped in seconds, so after
a failure, fix the cause and re-run the same command to resume.

Settings, from lowest to highest precedence:
  default config.yaml  <  --settings FILE  <  --set KEY=VALUE  <  dedicated flags

The code and default config are read from the workflow home: the folder
holding the workflow/ directory this script sits in (/opt/hemaclone inside
the image). Outputs, the sample sheet and logs belong to the project folder,
the current working directory.

Normally launched through Docker by hemaclone.sh; it can also run directly on a
machine where Snakemake and the tools are installed:
  python /path/to/workflow/scripts/run_batches.py --cores 16 --mem-gb 64 --batch-size 1
"""

import argparse
import copy
import csv
import datetime
import os
import shutil
import subprocess
import sys
import time

import yaml

WORKFLOW_HOME = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CONFIG = os.path.join(WORKFLOW_HOME, "config", "config.yaml")
PROFILE = os.path.join(WORKFLOW_HOME, "config", "profile")
SNAKEFILE = os.path.join(WORKFLOW_HOME, "workflow", "Snakefile")
VERSION_FILE = os.path.join(WORKFLOW_HOME, "VERSION")

# Project folder (current working directory)
BATCH_DIR = "config/batches"
RUN_LOG_DIR = "logs/batches"

# Snakemake options that the flags of this script control; passing them after
# "--" would bypass the checks below.
RESERVED_SNAKEMAKE_OPTIONS = ["--cores", "-c", "--jobs", "-j", "--resources", "--res",
                              "--config", "-C", "--configfile", "--configfiles",
                              "--snakefile", "-s", "--directory", "-d"]

# Keys that have a dedicated flag; --set on them is refused so that the
# launcher (which mounts folders and sets container limits) always sees them.
KEYS_WITH_FLAGS = {
    "machine.cores": "--cores",
    "machine.mem_gb": "--mem-gb",
    "batch_size": "--batch-size",
    "samples": "--sample-sheet",
    "reference.dir": "--resources-dir",
    "tmpdir": "--tmpdir",
    "fastq_dir": "--fastq-dir",
    "data_type": "--data-type",
}

SAMPLESHEET_SCRIPT = os.path.join(WORKFLOW_HOME, "workflow", "scripts", "make_samplesheet.py")


def parse_args(argv):
    if "--" in argv:
        split = argv.index("--")
        own_args, passthrough = argv[:split], argv[split + 1:]
    else:
        own_args, passthrough = argv, []

    parser = argparse.ArgumentParser(
        prog="run",
        description="Run HemaClone in batches. Every flag is optional; "
                    "defaults come from config/config.yaml.",
        epilog='Arguments after "--" are passed to snakemake unchanged '
               "(except cores, resources and config, which have their own flags).")
    parser.add_argument("--cores", type=int, default=None,
                        help="maximum cores for the run (default: machine.cores)")
    parser.add_argument("--mem-gb", type=float, default=None,
                        help="maximum RAM in GB for the run (default: machine.mem_gb)")
    parser.add_argument("--batch-size", type=int, default=None,
                        help="samples per batch (default: batch_size)")
    parser.add_argument("--only", default=None,
                        help="comma-separated sample names to run (default: every sample in the sheet)")
    parser.add_argument("--sample-sheet", default=None,
                        help="sample sheet to use (default: samples); built from the FASTQ folder "
                             "when it does not exist")
    parser.add_argument("--fastq-dir", default=None,
                        help="folder with the FASTQ files (default: fastq_dir)")
    parser.add_argument("--rebuild-sheet", action="store_true",
                        help="rebuild the sample sheet from the FASTQ folder, replacing the existing one")
    parser.add_argument("--allow-unparsed-headers", action="store_true",
                        help="when building the sheet, accept FASTQs whose first header is not "
                             "Casava 1.8+ (read group falls back to <sample>.<unit>)")
    parser.add_argument("--data-type", choices=["wes", "wgs"], default=None,
                        help="sequencing design (default: data_type); wgs needs "
                             "--set ch_filter.min_dp=... and --set ch_filter.min_vaf=...")
    parser.add_argument("--resources-dir", default=None,
                        help="folder with the reference files (default: reference.dir)")
    parser.add_argument("--tmpdir", default=None,
                        help="scratch folder for temporary files (default: tmpdir)")
    parser.add_argument("--settings", default=None, metavar="FILE",
                        help="YAML file with site-specific values, same keys as config/config.yaml")
    parser.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                        help="change any config value for this run, e.g. ch_filter.min_vaf=0.01 "
                             "or job_sizes.mutect2.mem_mb=8000; repeatable")
    parser.add_argument("--dry-run", action="store_true",
                        help="show the settings, batches and planned jobs; run nothing")
    parser.add_argument("--cohort-only", action="store_true",
                        help="skip the batches and only build the cohort QC report")
    args = parser.parse_args(own_args)

    reserved = [o for o in passthrough if o.split("=")[0] in RESERVED_SNAKEMAKE_OPTIONS]
    if reserved:
        sys.exit(f"ERROR: {reserved} cannot be passed after '--'. Use --cores, --mem-gb, "
                 "--set or --sample-sheet instead.")
    return args, passthrough


def fail(message):
    sys.exit(f"ERROR: {message}")


def preflight():
    problems = []
    for path in [CONFIG, os.path.join(PROFILE, "config.yaml"), SNAKEFILE]:
        if not os.path.isfile(path):
            problems.append(f"workflow file not found: {path}")
    if shutil.which("snakemake") is None:
        problems.append("snakemake not found on PATH")
    if problems:
        fail("\n  " + "\n  ".join(problems))


def workflow_version():
    """Version, code location and image tag, recorded with every run."""
    version = "unknown"
    if os.path.isfile(VERSION_FILE):
        with open(VERSION_FILE) as handle:
            version = handle.read().strip()
    info = {"version": version, "code": WORKFLOW_HOME,
            "image": os.environ.get("HEMACLONE_IMAGE", "none (not run through hemaclone.sh)"),
            "mode": os.environ.get("HEMACLONE_MODE", "direct")}
    if shutil.which("git"):
        commit = subprocess.run(["git", "-C", WORKFLOW_HOME, "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True)
        if commit.returncode == 0:
            dirty = subprocess.run(["git", "-C", WORKFLOW_HOME, "status", "--porcelain"],
                                   capture_output=True, text=True).stdout.strip()
            info["git_commit"] = commit.stdout.strip() + (" (with uncommitted changes)" if dirty else "")
    return info


# -----------------------------------------------------------------------------
# Overrides: --settings file and --set KEY=VALUE
# -----------------------------------------------------------------------------
def is_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def check_key(path, value, base_config, source):
    """Stop if the dotted key does not exist in config.yaml or the value has the wrong type."""
    key = ".".join(path)
    node = base_config
    for depth, part in enumerate(path):
        if not isinstance(node, dict) or part not in node:
            known = ", ".join(sorted(node)) if isinstance(node, dict) else "(none)"
            fail(f"{source} {key}: '{'.'.join(path[:depth + 1])}' is not a config key. "
                 f"Keys at this level: {known}")
        node = node[part]
    if isinstance(node, dict):
        fail(f"{source} {key}: this is a section, set one of its keys: {', '.join(sorted(node))}")
    if is_number(node) and not is_number(value):
        fail(f"{source} {key}: expected a number, got '{value}'")


def set_nested(target, path, value):
    for part in path[:-1]:
        target = target.setdefault(part, {})
    target[path[-1]] = value


def leaves(mapping, prefix=()):
    for key, value in mapping.items():
        if isinstance(value, dict):
            yield from leaves(value, prefix + (key,))
        else:
            yield prefix + (key,), value


def load_settings_file(path, base_config):
    if not os.path.isfile(path):
        fail(f"settings file not found: {path}")
    with open(path) as handle:
        settings = yaml.safe_load(handle) or {}
    if not isinstance(settings, dict):
        fail(f"settings file {path} must contain YAML keys like config/config.yaml")
    for key_path, value in leaves(settings):
        check_key(list(key_path), value, base_config, f"--settings {path}:")
    return settings


def parse_set_overrides(items, base_config):
    overrides = {}
    for item in items:
        if "=" not in item:
            fail(f"--set expects KEY=VALUE, got '{item}'")
        key, raw_value = item.split("=", 1)
        if key in KEYS_WITH_FLAGS:
            fail(f"use {KEYS_WITH_FLAGS[key]} instead of --set {key}")
        value = yaml.safe_load(raw_value)
        check_key(key.split("."), value, base_config, "--set")
        set_nested(overrides, key.split("."), value)
    return overrides


def deep_update(base, updates):
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = copy.deepcopy(value)
    return base


# -----------------------------------------------------------------------------
# Machine limits
# -----------------------------------------------------------------------------
def host_memory_gb():
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return None


def check_machine(cores, mem_gb):
    if not isinstance(cores, int) or cores < 1:
        fail(f"cores must be a whole number >= 1, got {cores!r}")
    if not is_number(mem_gb) or mem_gb <= 0:
        fail(f"mem_gb must be a positive number, got {mem_gb!r}")
    host_cores = os.cpu_count()
    if host_cores and cores > host_cores:
        fail(f"--cores {cores} is more than this machine has ({host_cores})")
    total_gb = host_memory_gb()
    if total_gb and mem_gb > total_gb:
        fail(f"--mem-gb {mem_gb:g} is more than this machine has ({total_gb:.0f} GB)")


# -----------------------------------------------------------------------------
# Sample sheet
# -----------------------------------------------------------------------------
def read_sheet(path):
    if not os.path.isfile(path):
        fail(f"sample sheet not found: {path}")
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not rows:
        fail(f"sample sheet is empty: {path}")
    samples = []
    for row in rows:
        if row["sample"] not in samples:
            samples.append(row["sample"])
    return fieldnames, rows, samples


def write_sheet(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_sheet(sheet, fastq_dir, rebuild, allow_unparsed):
    """Build the sample sheet from the FASTQ folder; stops the run if it cannot."""
    if not os.path.isdir(fastq_dir):
        fail(f"sample sheet {sheet} does not exist and FASTQ folder {fastq_dir} was not found "
             "(set it with --fastq-dir)")
    reason = "rebuilding it (--rebuild-sheet)" if rebuild else "it does not exist yet"
    print(f"--- sample sheet {sheet}: {reason}; building it from {fastq_dir}/", flush=True)
    command = [sys.executable, SAMPLESHEET_SCRIPT, "--fastq-dir", fastq_dir, "--out", sheet, "--force"]
    if allow_unparsed:
        command.append("--allow-unparsed-headers")
    if subprocess.run(command).returncode != 0:
        fail(f"could not build the sample sheet from {fastq_dir}/ (see the errors above)")


def warn_unlisted_fastq(rows, fastq_dir):
    """Warn about FASTQ, BAM or CRAM files in the folder that the sample sheet does not list."""
    if not os.path.isdir(fastq_dir):
        return
    listed = {os.path.realpath(row[col]) for row in rows for col in ("fq1", "fq2", "bam")
              if row.get(col)}
    unlisted = sorted(
        name for name in os.listdir(fastq_dir)
        if name.endswith((".fastq.gz", ".fq.gz", ".bam", ".cram")) and not name.startswith("Undetermined")
        and os.path.realpath(os.path.join(fastq_dir, name)) not in listed
    )
    if unlisted:
        print(f"WARNING: {len(unlisted)} input file(s) in {fastq_dir}/ are not in the sample sheet "
              "and will not be processed:", flush=True)
        for name in unlisted:
            print(f"    {name}", flush=True)
        print("  If they are new samples, add --rebuild-sheet (or rebuild with: hemaclone.sh samplesheet --force).",
              flush=True)


def resolved_reference(reference):
    folder = reference["dir"]
    return {key: (value if key == "dir" or os.path.isabs(value) else os.path.join(folder, value))
            for key, value in reference.items()}


def resolved_pon(config):
    """PoN path with the same lookup as the reference files, or None when turned off."""
    vcf = config["pon"].get("vcf")
    if not vcf:
        return None
    return vcf if os.path.isabs(vcf) else os.path.join(config["reference"]["dir"], vcf)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    args, passthrough = parse_args(sys.argv[1:])
    preflight()

    with open(CONFIG) as handle:
        base_config = yaml.safe_load(handle)

    # All overrides merged into one file that Snakemake reads after config.yaml
    overrides = {}
    if args.settings:
        deep_update(overrides, load_settings_file(args.settings, base_config))
    deep_update(overrides, parse_set_overrides(args.set, base_config))
    if args.resources_dir:
        set_nested(overrides, ["reference", "dir"], args.resources_dir)
    if args.tmpdir:
        set_nested(overrides, ["tmpdir"], args.tmpdir)
    if args.fastq_dir:
        set_nested(overrides, ["fastq_dir"], args.fastq_dir)
    if args.data_type:
        set_nested(overrides, ["data_type"], args.data_type)
    config = deep_update(copy.deepcopy(base_config), overrides)

    # The CH depth and VAF floors are tuned for WES depth: a WGS run must choose its own
    if config["data_type"] == "wgs":
        explicit = overrides.get("ch_filter", {})
        unset = [f"ch_filter.{key}" for key in ("min_dp", "min_vaf") if key not in explicit]
        if unset:
            fail(f"a WGS run must set {' and '.join(unset)} at launch (the defaults, DP >= "
                 f"{base_config['ch_filter']['min_dp']} and VAF >= {base_config['ch_filter']['min_vaf']}, "
                 "are tuned for WES depth), e.g. --set ch_filter.min_dp=20 --set ch_filter.min_vaf=0.10 "
                 "for 30x; see the README for the expected sensitivity")

    machine = config["machine"]
    cores = args.cores if args.cores is not None else machine["cores"]
    mem_gb = args.mem_gb if args.mem_gb is not None else machine["mem_gb"]
    check_machine(cores, mem_gb)
    budget_mb = int(mem_gb * 1024 * machine["schedule_fraction"])

    if not os.path.isdir(config["reference"]["dir"]):
        fail(f"resources folder not found: {config['reference']['dir']} (set it with --resources-dir)")

    # Sample sheet: built from the FASTQ folder when missing (or with --rebuild-sheet);
    # an existing sheet is used as it is, with a warning for FASTQs it does not list
    full_sheet = args.sample_sheet or config["samples"]
    fastq_dir = config["fastq_dir"]
    if args.rebuild_sheet or not os.path.isfile(full_sheet):
        build_sheet(full_sheet, fastq_dir, args.rebuild_sheet, args.allow_unparsed_headers)
    elif args.allow_unparsed_headers:
        print("NOTE: --allow-unparsed-headers only applies when the sheet is built "
              "(missing sheet or --rebuild-sheet).", flush=True)
    fieldnames, rows, all_samples = read_sheet(full_sheet)
    warn_unlisted_fastq(rows, fastq_dir)

    batch_size = args.batch_size if args.batch_size is not None else config["batch_size"]
    if not isinstance(batch_size, int) or batch_size < 1:
        fail(f"batch size must be a whole number >= 1, got {batch_size!r}")

    samples = all_samples
    if args.only:
        requested = [s.strip() for s in args.only.split(",") if s.strip()]
        unknown = [s for s in requested if s not in all_samples]
        if unknown:
            fail(f"samples not in {full_sheet}: {unknown}")
        samples = [s for s in all_samples if s in requested]
    batches = [samples[i:i + batch_size] for i in range(0, len(samples), batch_size)]

    # ---- record of this run's settings ----------------------------------------
    os.makedirs(BATCH_DIR, exist_ok=True)
    os.makedirs(RUN_LOG_DIR, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    run_prefix = os.path.join(RUN_LOG_DIR, f"run_{stamp}")

    override_file = None
    if overrides:
        override_file = f"{run_prefix}_overrides.yaml"
        with open(override_file, "w") as handle:
            yaml.safe_dump(overrides, handle, sort_keys=False)

    version = workflow_version()
    settings = {"workflow": version}
    settings.update(copy.deepcopy(config))
    settings["samples"] = full_sheet
    settings["machine"] = {"cores": cores, "mem_gb": mem_gb, "schedule_mb": budget_mb}
    settings["batch_size"] = batch_size
    settings["reference_resolved"] = resolved_reference(config["reference"])
    settings["pon_resolved"] = resolved_pon(config)
    settings["run"] = {"command": " ".join(sys.argv), "only": args.only,
                       "dry_run": args.dry_run, "snakemake_extra": passthrough}
    with open(f"{run_prefix}_settings.yaml", "w") as handle:
        yaml.safe_dump(settings, handle, sort_keys=False)

    pon = config["pon"]
    print("Settings for this run")
    print(f"  workflow     : {version['version']}  [{version['mode']}]  image: {version['image']}")
    print(f"  code from    : {WORKFLOW_HOME}" + (f"  (git {version['git_commit']})" if "git_commit" in version else ""))
    print(f"  sample sheet : {full_sheet}  (FASTQ folder: {fastq_dir})")
    print(f"  data type    : {config['data_type']}")
    print(f"  samples      : {len(samples)}   batch size: {batch_size}   batches: {len(batches)}")
    print(f"  cores        : {cores}")
    print(f"  RAM          : {mem_gb:g} GB limit, jobs scheduled into {budget_mb / 1024:.0f} GB")
    print(f"  resources    : {config['reference']['dir']}")
    print(f"  tmpdir       : {config['tmpdir']}")
    print(f"  CH filters   : DP >= {config['ch_filter']['min_dp']}, VAF >= {config['ch_filter']['min_vaf']}")
    pon_path = resolved_pon(config)
    print(f"  PoN          : {pon_path + ' (' + pon['mode'] + ')' if pon_path else 'none'}")
    print(f"  settings file: {args.settings or 'none'}")
    print(f"  --set        : {', '.join(args.set) if args.set else 'none'}")
    print(f"  full record  : {run_prefix}_settings.yaml")
    for number, batch in enumerate(batches, start=1):
        print(f"  batch {number:03d}: {', '.join(batch)}")

    def snakemake_command(targets, sheet, extra):
        command = ["snakemake", *targets, "--snakefile", SNAKEFILE, "--directory", os.getcwd(),
                   "--workflow-profile", PROFILE,
                   "--cores", str(cores), "--resources", f"mem_mb={budget_mb}"]
        command += extra
        if override_file:
            command += ["--configfile", override_file]
        command += ["--config", f"samples={sheet}"]
        command += passthrough
        return command

    cohort_extra = ["--allowed-rules", "multiqc", "qc_flags", "chip_cohort", "cohort_qc"]

    # ---- preflight once, on the whole sheet ----------------------------------------
    # The checks in common.smk run whenever Snakemake reads the workflow. A failure
    # there (a missing resource, a bad setting) concerns the run, not one sample:
    # stop before the first batch instead of marking every sample as failed.
    print("\n--- Preflight checks on the full sample sheet", flush=True)
    check = subprocess.run(snakemake_command([], full_sheet, ["--list-rules"]),
                           capture_output=True, text=True)
    if check.returncode != 0:
        print(check.stderr, file=sys.stderr)
        fail("preflight checks failed, so no batch was started. "
             "Fix the problem above and re-run the same command.")
    print("Preflight checks passed.", flush=True)

    # ---- cohort QC only ----------------------------------------------------------
    if args.cohort_only:
        extra = cohort_extra + (["--dry-run"] if args.dry_run else [])
        sys.exit(subprocess.run(snakemake_command(["cohort_qc"], full_sheet, extra)).returncode)

    run_log = f"{run_prefix}.tsv"
    with open(run_log, "w") as handle:
        handle.write("batch\tsamples\tstatus\tstarted\tminutes\n")

    def sheet_rows_of(sample_names):
        return [row for row in rows if row["sample"] in sample_names]

    coverage_metrics = "wgs_metrics" if config["data_type"] == "wgs" else "hs_metrics"

    def sample_done(sample):
        """A sample is finished when its two MAFs and its coverage metrics exist."""
        return all(os.path.isfile(path) for path in [
            f"maf/{sample}_all.funcotated.maf",
            f"maf/{sample}_ch_final.funcotated.maf",
            f"qc/bam/{sample}.{coverage_metrics}.txt",
        ])

    failed = {}     # sample -> reason, over the whole run

    # ---- run batches one after another; a failed sample does not stop the run ------
    for number, batch in enumerate(batches, start=1):
        started = datetime.datetime.now()
        t0 = time.time()
        print(f"\n=== batch {number:03d}/{len(batches):03d}  {', '.join(batch)}  "
              f"started {started:%Y-%m-%d %H:%M} ===", flush=True)
        batch_failed = {}

        # Step 1: integrity of the input files (FASTQ pairs, BAM/CRAM) of the samples not finished yet
        unfinished = [s for s in batch if not sample_done(s)]
        if unfinished and not args.dry_run:
            check_sheet = os.path.join(BATCH_DIR, f"batch_{number:03d}_check.tsv")
            check_units = sheet_rows_of(unfinished)
            write_sheet(check_sheet, fieldnames, check_units)
            checks = {}     # target -> (sample, label, log)
            for u in check_units:
                kind = "bam_check" if u.get("bam") else "fastq_check"
                name = f"{u['sample']}.{u['unit']}"
                label = "BAM/CRAM check" if u.get("bam") else "FASTQ check"
                checks[f"qc/{kind}/{name}.ok"] = (u["sample"], label, f"logs/{kind}/{name}.log")
            print(f"--- Input check: {len(checks)} unit(s) of {len(unfinished)} unfinished sample(s)", flush=True)
            subprocess.run(snakemake_command(list(checks), check_sheet, ["--quiet", "rules"]))
            for target, (sample, label, log_path) in checks.items():
                if not os.path.isfile(target):
                    batch_failed[sample] = f"{label} failed: {log_path}"

        # Step 2: the workflow, for every sample whose FASTQs passed
        to_run = [s for s in batch if s not in batch_failed]
        if to_run:
            sheet = os.path.join(BATCH_DIR, f"batch_{number:03d}.tsv")
            write_sheet(sheet, fieldnames, sheet_rows_of(to_run))
            extra = ["--dry-run"] if args.dry_run else []
            returncode = subprocess.run(snakemake_command(["per_sample"], sheet, extra)).returncode
            for sample in to_run:
                if args.dry_run:
                    if returncode != 0:
                        batch_failed[sample] = f"dry run failed (exit {returncode})"
                elif not sample_done(sample):
                    batch_failed[sample] = ("workflow step failed: "
                                            f"logs/<step>/{sample}*.log, .snakemake/log/")

        minutes = (time.time() - t0) / 60
        if batch_failed:
            status = "failed: " + ", ".join(sorted(batch_failed))
        else:
            status = "ok"
        failed.update(batch_failed)

        with open(run_log, "a") as handle:
            handle.write(f"{number:03d}\t{','.join(batch)}\t{status}\t{started:%Y-%m-%d %H:%M}\t{minutes:.1f}\n")
        print(f"=== batch {number:03d}: {status}, {minutes:.1f} min ===", flush=True)
        for sample, reason in sorted(batch_failed.items()):
            print(f"    {sample}: {reason}", flush=True)

    # ---- summary ----------------------------------------------------------------------
    if failed:
        failed_file = os.path.join(RUN_LOG_DIR, f"failed_{stamp}.tsv")
        with open(failed_file, "w") as handle:
            handle.write("sample\treason\n")
            for sample, reason in sorted(failed.items()):
                handle.write(f"{sample}\t{reason}\n")
        print(f"\n{len(failed)} of {len(samples)} sample(s) failed; the others are finished:", file=sys.stderr)
        for sample, reason in sorted(failed.items()):
            print(f"  {sample}: {reason}", file=sys.stderr)
        print(f"\n  List: {failed_file}\n"
              f"  Fix the cause and re-run the same command: finished samples are skipped.\n"
              f"  Cohort QC report not built: it runs once every sample has finished.\n"
              f"  Batch summary: {run_log}", file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print("\nDry run: FASTQ checks and the cohort QC step run only in a real run.")
        return

    # ---- cohort QC report over the whole sheet -------------------------------------
    if args.only:
        print("\nCohort QC report not rebuilt because --only was used. "
              "Once every sample is done, run again with --cohort-only.")
        print(f"Batch summary: {run_log}")
        return

    print("\n=== cohort QC report ===", flush=True)
    returncode = subprocess.run(snakemake_command(["cohort_qc"], full_sheet, cohort_extra)).returncode
    with open(run_log, "a") as handle:
        handle.write(f"cohort\tall\t{'ok' if returncode == 0 else 'failed'}\t-\t-\n")
    print(f"Batch summary: {run_log}")
    sys.exit(returncode)


if __name__ == "__main__":
    main()
