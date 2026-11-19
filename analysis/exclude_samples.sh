#!/usr/bin/env bash
# =============================================================================
# Remove samples from the HemaClone cohort tables read by ch_analysis.Rmd.
#
# The original chip/ folder is copied once to chip_unfiltered/ and never
# modified; the four cohort tables in chip/ are rewritten from those originals,
# so the script can be rerun (e.g. with a different list) without stacking
# filters. Per-sample files of excluded samples are removed from chip/samples/
# (the originals stay in chip_unfiltered/samples/).
#
# Usage: bash exclude_samples.sh --exclude IDS_FILE PROJECT_DIR
#   IDS_FILE     one sample ID per line; blank lines and lines starting with #
#                are ignored. Keep it with the data, not in the repository.
#   PROJECT_DIR  folder holding chip/ (the workflow project folder)
# =============================================================================
set -euo pipefail

usage() {
  echo "Usage: bash exclude_samples.sh --exclude IDS_FILE PROJECT_DIR" >&2
  exit 1
}

IDS_FILE=""
PROJECT_DIR=""
while (( $# > 0 )); do
  case "$1" in
    --exclude) [[ $# -ge 2 ]] || usage; IDS_FILE="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) [[ -z "$PROJECT_DIR" ]] || usage; PROJECT_DIR="$1"; shift ;;
  esac
done
[[ -n "$IDS_FILE" && -n "$PROJECT_DIR" ]] || usage
[[ -f "$IDS_FILE" ]] || { echo "ERROR: IDs file not found: $IDS_FILE" >&2; exit 1; }
[[ -d "$PROJECT_DIR" ]] || { echo "ERROR: project folder not found: $PROJECT_DIR" >&2; exit 1; }

# IDs: one per line, trimmed, without blanks and comments
EXCLUDE=()
while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%%#*}"
  line="$(echo "$line" | tr -d '[:space:]')"
  [[ -n "$line" ]] && EXCLUDE+=("$line")
done < "$IDS_FILE"
if (( ${#EXCLUDE[@]} == 0 )); then
  echo "ERROR: no sample IDs in $IDS_FILE" >&2
  exit 1
fi
echo "Samples to exclude (${#EXCLUDE[@]}): ${EXCLUDE[*]}"

TABLES=(chip_samples.tsv chip_calls.maf chip_calls.oncoplot.maf chip_candidates.tsv)
BACKUP="chip_unfiltered"

cd "$PROJECT_DIR"

# ---- preflight ----------------------------------------------------------------
for t in "${TABLES[@]}"; do
  if [[ ! -f "chip/$t" && ! -f "$BACKUP/$t" ]]; then
    echo "ERROR: chip/$t not found in $PROJECT_DIR" >&2
    exit 1
  fi
done

if [[ ! -d "$BACKUP" ]]; then
  cp -a chip "$BACKUP"
  echo "Originals copied to $PROJECT_DIR/$BACKUP/"
else
  echo "Using the originals already in $PROJECT_DIR/$BACKUP/"
fi

# Every sample to exclude must exist in the original sample table
missing=()
for s in "${EXCLUDE[@]}"; do
  if ! cut -f1 "$BACKUP/chip_samples.tsv" | grep -qxF -- "$s"; then
    missing+=("$s")
  fi
done
if (( ${#missing[@]} > 0 )); then
  echo "ERROR: samples not found in $BACKUP/chip_samples.tsv: ${missing[*]}" >&2
  exit 1
fi

EXCLUDE_LIST=$(IFS=,; echo "${EXCLUDE[*]}")

# ---- filter: keeps comment lines and the header, drops rows of excluded samples
filter_table() {
  local in="$1" out="$2"
  awk -F'\t' -v ex="$EXCLUDE_LIST" '
    BEGIN { n = split(ex, a, ","); for (i = 1; i <= n; i++) drop[a[i]] = 1 }
    !header && /^#/ { print; next }
    !header {
      header = 1
      for (i = 1; i <= NF; i++) if ($i == "Tumor_Sample_Barcode") col = i
      if (!col) { bad = 1; exit }
      print; next
    }
    $col in drop { removed++; next }
    { kept++; print }
    END {
      if (bad) { print "no Tumor_Sample_Barcode column" > "/dev/stderr"; exit 2 }
      printf "%-28s kept %6d rows, removed %6d\n", FILENAME, kept, removed > "/dev/stderr"
    }' "$in" > "$out.tmp"
  mv "$out.tmp" "$out"
}

for t in "${TABLES[@]}"; do
  filter_table "$BACKUP/$t" "chip/$t"
done

# ---- per-sample files ---------------------------------------------------------
for s in "${EXCLUDE[@]}"; do
  rm -f "chip/samples/$s".*
done

# ---- verification -------------------------------------------------------------
for t in "${TABLES[@]}"; do
  for s in "${EXCLUDE[@]}"; do
    if grep -qF -- "$s" "chip/$t"; then
      echo "ERROR: $s still present in chip/$t" >&2
      exit 1
    fi
  done
done

n_before=$(( $(wc -l < "$BACKUP/chip_samples.tsv") - 1 ))
n_after=$(( $(wc -l < chip/chip_samples.tsv) - 1 ))
echo "Samples: $n_before -> $n_after (excluded ${#EXCLUDE[@]})"
if (( n_before - n_after != ${#EXCLUDE[@]} )); then
  echo "ERROR: expected to remove ${#EXCLUDE[@]} samples, removed $(( n_before - n_after ))" >&2
  exit 1
fi
echo "Done. chip/ now holds the filtered tables; originals in $BACKUP/"
