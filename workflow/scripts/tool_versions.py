# Versions of the tools the run used, for the methods report (Snakemake script).
#
# The pinned versions come from the environment files the image is built from
# (workflow/docker/*.yaml); the installed versions from the conda metadata of
# each environment: environment.yaml is installed into the base environment
# (the one running this script), the others into envs/<name>. Outside the
# image the environments may not exist: their versions are written as
# "not found", so the table shows what was checked.

import glob
import json
import os
import sys

import yaml

# Environment file -> environment folder name; None = the base environment
ENVIRONMENTS = [("environment.yaml", None), ("somalier.yaml", "somalier"), ("report.yaml", "report")]


def installed_versions(prefix):
    """Package name -> version, from <prefix>/conda-meta/*.json."""
    versions = {}
    for path in glob.glob(os.path.join(prefix, "conda-meta", "*.json")):
        with open(path) as handle:
            record = json.load(handle)
        if "name" in record and "version" in record:
            versions[record["name"]] = record["version"]
    return versions


def main(snakemake, log):
    docker_dir = snakemake.params.docker_dir
    base_prefix = sys.prefix
    rows = []
    for file_name, env_name in ENVIRONMENTS:
        path = os.path.join(docker_dir, file_name)
        if not os.path.isfile(path):
            print(f"ERROR: environment file not found: {path}", file=log, flush=True)
            sys.exit(f"environment file not found: {path}")
        with open(path) as handle:
            spec = yaml.safe_load(handle)
        prefix = base_prefix if env_name is None else os.path.join(base_prefix, "envs", env_name)
        installed = installed_versions(prefix)
        label = "base" if env_name is None else env_name
        print(f"{file_name}: environment {label} at {prefix}, {len(installed)} packages installed", file=log)
        for dependency in spec["dependencies"]:
            name, _, pinned = str(dependency).partition("=")
            rows.append((label, name, pinned or "-", installed.get(name, "not found")))

    with open(snakemake.output[0], "w") as handle:
        handle.write("environment\tpackage\tpinned\tinstalled\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    # A pin such as python=3.12 matches any 3.12.x
    differing = [r for r in rows if r[3] != "not found" and r[3] != r[2]
                 and not r[3].startswith(r[2] + ".")]
    print(f"packages: {len(rows)} | installed version differs from the pin: {len(differing)}", file=log)
    for row in differing:
        print(f"  {row[0]} {row[1]}: pinned {row[2]}, installed {row[3]}", file=log)


with open(snakemake.log[0], "w") as log_handle:
    main(snakemake, log_handle)
