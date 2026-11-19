#!/usr/bin/env bash
# =============================================================================
# Check that two FASTQ mates are intact and in step, in one pass over both:
#   - both files decompress without error
#   - both hold only complete 4-line records
#   - both hold the same number of reads
#   - read i has the same name in both files (comment and /1 /2 removed)
# Writes <output> only when every check passes; otherwise exits 1 with the
# first problem found.
#
# Usage: check_fastq_pair.sh <R1.fastq.gz> <R2.fastq.gz> <output>
# =============================================================================

set -euo pipefail

r1="$1"
r2="$2"
out="$3"

work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# Read names, one per record, then a final "##END <lines>" line.
# The zcat exit status is saved to a file, since a process substitution's
# status is not seen by the main shell.
names() {
    local fastq="$1" status_file="$2"
    {
        if zcat "$fastq"; then status=0; else status=$?; fi
        echo "$status" > "$status_file"
    } | awk 'NR % 4 == 1 { name = $1; sub(/\/[12]$/, "", name); print name }
             END { print "##END " NR }'
}

echo "Checking $r1 and $r2"

compare_status=0
paste <(names "$r1" "$work/r1.status") <(names "$r2" "$work/r2.status") \
    | awk -F'\t' '
        $1 ~ /^##END / && $2 ~ /^##END / {
            lines1 = substr($1, 7) + 0; lines2 = substr($2, 7) + 0
            if (lines1 % 4 != 0) { print "R1 ends with an incomplete record (" lines1 " lines)"; exit 1 }
            if (lines2 % 4 != 0) { print "R2 ends with an incomplete record (" lines2 " lines)"; exit 1 }
            print "reads\t" lines1 / 4
            exit 0
        }
        $1 ~ /^##END / || $2 ~ /^##END / {
            short = ($1 ~ /^##END /) ? "R1" : "R2"
            print "read counts differ: " short " ends at read " NR - 1 " while the other file continues"
            exit 1
        }
        $1 != $2 {
            print "mates out of step at read " NR ": R1 \"" $1 "\" vs R2 \"" $2 "\""
            exit 1
        }' > "$work/result" || compare_status=$?

cat "$work/result"

# 141 = zcat stopped because the comparison had already failed (broken pipe),
# not a problem with that file
gzip_problem=""
for mate in r1 r2; do
    if [[ -f "$work/$mate.status" ]]; then
        status="$(cat "$work/$mate.status")"
        if [[ "$status" != "0" && "$status" != "141" ]]; then
            gzip_problem+="${mate^^} does not decompress cleanly (zcat exit $status). "
        fi
    fi
done

if [[ -n "$gzip_problem" ]]; then
    echo "$gzip_problem"
fi

if [[ $compare_status -ne 0 || -n "$gzip_problem" ]] || ! grep -q "^reads" "$work/result"; then
    echo "FAILED: $r1 / $r2 are damaged or out of step. Re-copy them from the source."
    exit 1
fi

{
    echo -e "r1\t$r1"
    echo -e "r2\t$r2"
    grep "^reads" "$work/result"
    echo -e "checked\t$(date '+%Y-%m-%d %H:%M:%S')"
} > "$out"
echo "OK"
