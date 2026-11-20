#!/bin/bash
# =============================================================================
# Checks run by CI on every pull request; also runnable on any machine with
# Docker, against an image built from the checkout being tested:
#
#   bash workflow/hemaclone.sh build
#   bash tests/ci_checks.sh hemaclone:<VERSION>
#
# 1. Syntax of every shell and Python script in the image
# 2. Synthetic project (tests/make_fixture.sh): tiny FASTQs, empty resources
# 3. Full launcher dry run, with the launcher copied out of the image
# 4. Dry run of the whole job graph, cohort steps included
# 5. Snakemake lint: fails on any lint but the per-rule container lint
#
# Nothing is aligned or called: these checks catch broken rules, settings,
# preflight checks and launcher logic, not wrong results.
# =============================================================================

set -euo pipefail

[[ $# -eq 1 ]] || { echo "usage: bash tests/ci_checks.sh IMAGE" >&2; exit 1; }
IMAGE="$1"

CORES=4
MEM_GB=14
# Jobs are scheduled into 0.9 of the memory limit (machine.schedule_fraction)
BUDGET_MB=$(( MEM_GB * 1024 * 9 / 10 ))
IMAGE_HOME="/opt/hemaclone"
TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: image $IMAGE not found" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
PROJECT="${WORK}/project"
RESOURCES="${WORK}/resources"
TMP_DIR="${WORK}/tmp"
mkdir -p "$PROJECT" "$TMP_DIR"

echo "== 1. Script syntax"
docker run --rm "$IMAGE" bash -c "
    set -e
    for f in ${IMAGE_HOME}/workflow/*.sh ${IMAGE_HOME}/workflow/scripts/*.sh; do bash -n \"\$f\"; done
    python -c 'import ast, sys
for path in sys.argv[1:]:
    with open(path) as handle:
        ast.parse(handle.read(), path)
print(f\"{len(sys.argv) - 1} Python scripts parsed\")' ${IMAGE_HOME}/workflow/scripts/*.py
"

echo "== 2. Synthetic project"
bash "${TESTS_DIR}/make_fixture.sh" "$PROJECT" "$RESOURCES"

echo "== 3. Launcher dry run"
cd "$PROJECT"
docker run --rm "$IMAGE" cat "${IMAGE_HOME}/workflow/hemaclone.sh" > hemaclone.sh
bash hemaclone.sh run \
    --image "$IMAGE" \
    --cores "$CORES" --mem-gb "$MEM_GB" --batch-size 1 \
    --resources-dir "$RESOURCES" \
    --tmpdir "$TMP_DIR" \
    --set job_sizes.bwa_map.mem_mb=12000 \
    --set job_sizes.markdup.mem_mb=12000 \
    --dry-run

OVERRIDES="$(ls logs/batches/run_*_overrides.yaml | tail -n 1)"
[[ -f "$OVERRIDES" ]] || { echo "ERROR: the dry run wrote no overrides file" >&2; exit 1; }

# snakemake inside the image, on the same project and settings as the dry run
snakemake_in_image() {
    docker run --rm \
        -u "$(id -u):$(id -g)" \
        -v "${PROJECT}:${PROJECT}" -v "${RESOURCES}:${RESOURCES}" -v "${TMP_DIR}:${TMP_DIR}" \
        -w "$PROJECT" \
        "$IMAGE" \
        snakemake "$@" \
            --snakefile "${IMAGE_HOME}/workflow/Snakefile" \
            --directory "$PROJECT" \
            --workflow-profile "${IMAGE_HOME}/config/profile" \
            --cores "$CORES" --resources "mem_mb=${BUDGET_MB}" \
            --configfile "$OVERRIDES" \
            --config samples=config/samples.tsv
}

echo "== 4. Whole job graph, cohort steps included"
snakemake_in_image all --dry-run --quiet rules

echo "== 5. Lint"
LINT_REPORT="${WORK}/lint.txt"
LINT_STATUS=0
snakemake_in_image --lint > "$LINT_REPORT" 2>&1 || LINT_STATUS=$?
grep -v '^\[preflight\]' "$LINT_REPORT" || true
# Snakemake exits 1 when it finds lints; a failure without any lint means lint did not run
if [[ "$LINT_STATUS" -ne 0 ]] && ! grep -q '^Lints for' "$LINT_REPORT"; then
    echo "ERROR: snakemake --lint failed (exit ${LINT_STATUS}) without reporting lints" >&2
    exit 1
fi
# The per-rule container lint does not apply: the whole workflow runs in one image
OTHER_LINTS="$(grep -E '^\s+\* ' "$LINT_REPORT" | grep -vc 'Specify a conda environment or container' || true)"
echo "Lints other than the per-rule container lint: ${OTHER_LINTS}"
if [[ "$OTHER_LINTS" -ne 0 ]]; then
    echo "ERROR: fix the lints above" >&2
    exit 1
fi

echo "All checks passed."
