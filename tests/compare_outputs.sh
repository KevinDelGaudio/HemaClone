#!/bin/bash
# =============================================================================
# Compare the results of two runs of the same samples, e.g. before and after a
# change that must leave every output identical:
#
#   bash tests/compare_outputs.sh OLD_PROJECT NEW_PROJECT
#
# Compared, as records (header lines starting with # skipped, since they hold
# dates, command lines and versions):
#   vcf/pass/*.vcf.gz, vcf/ch/*.ch_final.vcf.gz, vcf/gnomad/*.vcf.gz
#   maf/*.maf
#   chip/*.tsv, chip/*.maf, chip/samples/*.tsv
#   qc/qc_flags.tsv
# Not compared: logs, benchmarks, BAMs, QC reports with run dates.
#
# For each file: identical, or the number of records found only in the old run
# and only in the new one. Exits 1 if any file differs or exists in one run only.
# =============================================================================

set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: bash tests/compare_outputs.sh OLD_PROJECT NEW_PROJECT" >&2; exit 1; }
OLD="$1"
NEW="$2"
for dir in "$OLD" "$NEW"; do
    [[ -d "$dir" ]] || { echo "ERROR: folder not found: $dir" >&2; exit 1; }
done

PATTERNS=(
    "vcf/pass/*.vcf.gz"
    "vcf/ch/*.ch_final.vcf.gz"
    "vcf/gnomad/*.vcf.gz"
    "maf/*.maf"
    "chip/*.tsv"
    "chip/*.maf"
    "chip/samples/*.tsv"
    "qc/qc_flags.tsv"
)

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Records of a file, without header lines, sorted
records() {
    local path="$1"
    if [[ "$path" == *.gz ]]; then
        zcat "$path" | grep -v '^#' | LC_ALL=C sort || true
    else
        grep -v '^#' "$path" | LC_ALL=C sort || true
    fi
}

# Relative paths of the files matching the patterns in one project
list_files() {
    local project="$1" pattern path
    for pattern in "${PATTERNS[@]}"; do
        for path in "$project"/$pattern; do
            if [[ -f "$path" ]]; then
                echo "${path#"$project"/}"
            fi
        done
    done
}

list_files "$OLD" | sort > "$WORK/old_files"
list_files "$NEW" | sort > "$WORK/new_files"

n_old="$(wc -l < "$WORK/old_files")"
n_new="$(wc -l < "$WORK/new_files")"
echo "Files: ${n_old} in ${OLD}, ${n_new} in ${NEW}"
[[ "$n_old" -gt 0 ]] || { echo "ERROR: no result files in ${OLD}" >&2; exit 1; }

problems=0

while read -r file; do
    echo "ONLY IN OLD  ${file}"
    problems=$((problems + 1))
done < <(comm -23 "$WORK/old_files" "$WORK/new_files")

while read -r file; do
    echo "ONLY IN NEW  ${file}"
    problems=$((problems + 1))
done < <(comm -13 "$WORK/old_files" "$WORK/new_files")

identical=0
while read -r file; do
    records "$OLD/$file" > "$WORK/a"
    records "$NEW/$file" > "$WORK/b"
    rows="$(wc -l < "$WORK/a")"
    if cmp -s "$WORK/a" "$WORK/b"; then
        identical=$((identical + 1))
        printf 'identical    %s (%s records)\n' "$file" "$rows"
    else
        only_old="$(comm -23 "$WORK/a" "$WORK/b" | wc -l)"
        only_new="$(comm -13 "$WORK/a" "$WORK/b" | wc -l)"
        printf 'DIFFERENT    %s (%s records only in old, %s only in new)\n' "$file" "$only_old" "$only_new"
        problems=$((problems + 1))
    fi
done < <(comm -12 "$WORK/old_files" "$WORK/new_files")

echo
echo "Identical files: ${identical}; differences: ${problems}"
if [[ "$problems" -ne 0 ]]; then
    exit 1
fi
echo "Outputs identical."
