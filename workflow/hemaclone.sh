#!/bin/bash
# =============================================================================
# HemaClone launcher. Run it from the project folder (FASTQs in fastq/).
#
#   bash hemaclone.sh run [flags]              run everything: builds config/samples.tsv
#                                         from fastq/ when it does not exist yet
#   bash hemaclone.sh run --help               list every run flag
#   bash hemaclone.sh samplesheet              only build config/samples.tsv, to inspect it
#   bash hemaclone.sh shell                    shell inside the container
#   bash hemaclone.sh build                    build the image (from a code checkout)
#
# The workflow code and its defaults come from the image (--image, default:
# the image this launcher was copied from). --dev uses the code of the checkout
# this script sits in instead, with the tools of the same image.
#
# --cores and --mem-gb are hard limits: the container can never use more CPU
# or RAM than that. The project folder, the FASTQ folder, the resources folder
# and the tmpdir are mounted automatically; add other host folders with --mount PATH.
#
# Get the launcher for a new project from the image:
#   docker run --rm hemaclone:<version> cat /opt/hemaclone/workflow/hemaclone.sh > hemaclone.sh
# =============================================================================

set -euo pipefail

# Filled in when the image is built (see the Dockerfile)
DEFAULT_IMAGE=""

IMAGE_HOME="/opt/hemaclone"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHECKOUT="$(dirname "$SCRIPT_DIR")"
ROOT="$(pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }

# Help: the header comment, from the line after the first "# ===" rule to the next one
usage() {
    awk 'NR > 2 && /^# =+$/ { exit }
         NR > 2 { sub(/^# ?/, ""); print }' "${BASH_SOURCE[0]}"
    exit 0
}

# Value of a key in a simple YAML file: read_yaml FILE SECTION KEY
# (SECTION empty for a top-level key). Prints nothing if absent.
read_yaml() {
    local file="$1" section="$2" key="$3"
    [[ -n "$file" && -f "$file" ]] || return 0
    awk -v section="$section" -v key="$key" '
        function clean(v) {
            sub(/^[^:]*:[ \t]*/, "", v); sub(/[ \t]+#.*$/, "", v)
            gsub(/^["\047]|["\047]$/, "", v); return v
        }
        section == "" && $0 ~ "^" key ":"     { print clean($0); exit }
        section != "" && $0 ~ "^" section ":" { inside = 1; next }
        inside && /^[^ #]/                    { inside = 0 }
        inside && $1 == key ":"               { print clean($0); exit }
    ' "$file"
}

# Setting with precedence: flag value > settings file > default config
setting() {
    local flag_value="$1" section="$2" key="$3" value
    if [[ -n "$flag_value" ]]; then echo "$flag_value"; return; fi
    value="$(read_yaml "$SETTINGS" "$section" "$key")"
    [[ -n "$value" ]] || value="$(read_yaml "$DEFAULT_CONFIG" "$section" "$key")"
    echo "$value"
}

absolute() {
    if [[ "$1" == /* ]]; then echo "$1"; else echo "${ROOT}/$1"; fi
}

checkout_version() {
    [[ -f "${CHECKOUT}/VERSION" ]] || die "no VERSION file in ${CHECKOUT}: build and --dev need a code checkout"
    tr -d '[:space:]' < "${CHECKOUT}/VERSION"
}

resolve_image() {
    if [[ -n "$IMAGE_FLAG" ]]; then IMAGE="$IMAGE_FLAG"
    elif [[ -n "${HEMACLONE_IMAGE:-}" ]]; then IMAGE="$HEMACLONE_IMAGE"
    elif [[ -n "$DEFAULT_IMAGE" ]]; then IMAGE="$DEFAULT_IMAGE"
    elif [[ -f "${CHECKOUT}/VERSION" ]]; then IMAGE="hemaclone:$(checkout_version)"
    else die "no image given: use --image hemaclone:<version>"
    fi
    docker image inspect "$IMAGE" >/dev/null 2>&1 \
        || die "image $IMAGE not found. Build it from a code checkout: bash workflow/hemaclone.sh build"
}

# Where the code and default config come from: the image, or the checkout (--dev)
resolve_code() {
    if [[ $DEV -eq 1 ]]; then
        [[ -f "${CHECKOUT}/workflow/Snakefile" ]] || die "--dev: no workflow code in ${CHECKOUT}"
        MODE="dev"
        CODE_HOME="$CHECKOUT"
        DEFAULT_CONFIG="${CHECKOUT}/config/config.yaml"
    else
        MODE="image"
        CODE_HOME="$IMAGE_HOME"
        DEFAULT_CONFIG="$(mktemp)"
        trap 'rm -f "$DEFAULT_CONFIG"' EXIT
        docker run --rm "$IMAGE" cat "${IMAGE_HOME}/config/config.yaml" > "$DEFAULT_CONFIG" 2>/dev/null \
            || die "image $IMAGE contains no workflow code (built before the code was added?). Rebuild it with: bash workflow/hemaclone.sh build, or use --dev"
    fi
}

# FASTQ folder as the container must see it: its real location, so a fastq/
# that is a symlink to another drive is mounted too. Empty if it does not exist.
fastq_mount() {
    local fastq_dir
    fastq_dir="$(absolute "$(setting "$FASTQ_DIR" "" fastq_dir)")"
    [[ -d "$fastq_dir" ]] && readlink -f "$fastq_dir" || true
}

# Folders to mount: project, FASTQs, resources, tmpdir, code checkout (--dev), extras
build_mounts() {
    local resources_dir tmp_dir fastq_real m
    resources_dir="$(absolute "$(setting "$RESOURCES_DIR" reference dir)")"
    tmp_dir="$(absolute "$(setting "$TMP_DIR" "" tmpdir)")"
    fastq_real="$(fastq_mount)"

    [[ -d "$resources_dir" ]] || die "resources folder not found: $resources_dir (set it with --resources-dir)"
    mkdir -p "$tmp_dir" || die "cannot create tmpdir: $tmp_dir"

    MOUNTS=( "$ROOT" "$resources_dir" "$tmp_dir" "${EXTRA_MOUNTS[@]}" )
    [[ -n "$fastq_real" ]] && MOUNTS+=( "$fastq_real" )
    [[ "$MODE" == "dev" ]] && MOUNTS+=( "$CHECKOUT" )
    VOLUME_ARGS=()
    local seen=" "
    for m in "${MOUNTS[@]}"; do
        m="$(absolute "$m")"
        [[ -d "$m" ]] || die "folder to mount not found: $m"
        [[ "$seen" == *" $m "* ]] && continue
        seen+="$m "
        VOLUME_ARGS+=( -v "${m}:${m}" )
    done
}

# docker run with mounts, host user, and CPU/RAM hard limits
docker_run() {
    local cores="$1" mem_gb="$2"; shift 2
    local tty_args=()
    [[ -t 0 && -t 1 ]] && tty_args=( -it )

    docker run --rm --init "${tty_args[@]}" \
        --cpus "$cores" \
        --memory "${mem_gb}g" \
        --memory-swap "${mem_gb}g" \
        -u "$(id -u):$(id -g)" \
        -e HEMACLONE_IMAGE="$IMAGE" \
        -e HEMACLONE_MODE="$MODE" \
        "${VOLUME_ARGS[@]}" \
        -w "$ROOT" \
        "$IMAGE" \
        "$@"
}

[[ $# -ge 1 ]] || usage
COMMAND="$1"; shift

CORES=""; MEM_GB=""; RESOURCES_DIR=""; TMP_DIR=""; FASTQ_DIR=""; SETTINGS=""; IMAGE_FLAG=""; DEV=0
EXTRA_MOUNTS=(); RUN_ARGS=()

# Flags read by the launcher. --cores, --mem-gb, --mount, --image and --dev are
# consumed here; --resources-dir, --tmpdir, --fastq-dir and --settings are also passed on.
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores)         [[ $# -ge 2 ]] || die "--cores needs a value";          CORES="$2"; shift 2 ;;
        --mem-gb)        [[ $# -ge 2 ]] || die "--mem-gb needs a value";         MEM_GB="$2"; shift 2 ;;
        --mount)         [[ $# -ge 2 ]] || die "--mount needs a folder";         EXTRA_MOUNTS+=( "$2" ); shift 2 ;;
        --image)         [[ $# -ge 2 ]] || die "--image needs a name";           IMAGE_FLAG="$2"; shift 2 ;;
        --dev)           DEV=1; shift ;;
        --resources-dir) [[ $# -ge 2 ]] || die "--resources-dir needs a folder"; RESOURCES_DIR="$2"; RUN_ARGS+=( "$1" "$2" ); shift 2 ;;
        --tmpdir)        [[ $# -ge 2 ]] || die "--tmpdir needs a folder";        TMP_DIR="$2"; RUN_ARGS+=( "$1" "$2" ); shift 2 ;;
        --fastq-dir)     [[ $# -ge 2 ]] || die "--fastq-dir needs a folder";     FASTQ_DIR="$2"; RUN_ARGS+=( "$1" "$2" ); shift 2 ;;
        --settings)      [[ $# -ge 2 ]] || die "--settings needs a file";        SETTINGS="$2"; RUN_ARGS+=( "$1" "$2" ); shift 2 ;;
        --)              RUN_ARGS+=( "$@" ); break ;;
        *)               RUN_ARGS+=( "$1" ); shift ;;
    esac
done

[[ -z "$SETTINGS" || -f "$SETTINGS" ]] || die "settings file not found: $SETTINGS"

case "$COMMAND" in
    -h|--help|help)
        usage
        ;;

    build)
        VERSION="$(checkout_version)"
        TAG="${IMAGE_FLAG:-hemaclone:${VERSION}}"
        echo "Building ${TAG} from ${CHECKOUT}"
        docker build -f "${CHECKOUT}/workflow/docker/Dockerfile" -t "$TAG" "${RUN_ARGS[@]}" "$CHECKOUT"
        ;;

    samplesheet)
        resolve_image; resolve_code
        FASTQ_DIR="$(setting "$FASTQ_DIR" "" fastq_dir)"
        FASTQ_REAL="$(fastq_mount)"
        [[ -n "$FASTQ_REAL" ]] || die "FASTQ folder not found: $(absolute "$FASTQ_DIR") (set it with --fastq-dir)"
        VOLUME_ARGS=( -v "${ROOT}:${ROOT}" )
        [[ "$FASTQ_REAL" != "$ROOT" && "$FASTQ_REAL" != "$ROOT"/* ]] && VOLUME_ARGS+=( -v "${FASTQ_REAL}:${FASTQ_REAL}" )
        [[ "$MODE" == "dev" && "$CHECKOUT" != "$ROOT" ]] && VOLUME_ARGS+=( -v "${CHECKOUT}:${CHECKOUT}" )
        # make_samplesheet.py takes --fastq-dir, --out, --force, --allow-unparsed-headers
        docker_run 2 8 python "${CODE_HOME}/workflow/scripts/make_samplesheet.py" \
            --fastq-dir "$FASTQ_DIR" --out config/samples.tsv "${RUN_ARGS[@]}"
        ;;

    shell)
        resolve_image; resolve_code; build_mounts
        docker_run "$(setting "$CORES" machine cores)" "$(setting "$MEM_GB" machine mem_gb)" bash
        ;;

    run)
        resolve_image; resolve_code

        CORES="$(setting "$CORES" machine cores)"
        MEM_GB="$(setting "$MEM_GB" machine mem_gb)"
        [[ "$CORES"  =~ ^[0-9]+$ && "$CORES"  -ge 1 ]] || die "--cores must be a whole number >= 1 (got '$CORES')"
        [[ "$MEM_GB" =~ ^[0-9]+$ && "$MEM_GB" -ge 1 ]] || die "--mem-gb must be a whole number of GB >= 1 (got '$MEM_GB')"

        HOST_CORES="$(nproc)"
        HOST_MEM_GB="$(awk '/^MemTotal:/ { printf "%d", $2 / 1024 / 1024 }' /proc/meminfo)"
        [[ "$CORES"  -le "$HOST_CORES" ]]  || die "--cores $CORES is more than this machine has ($HOST_CORES)"
        [[ "$MEM_GB" -le "$HOST_MEM_GB" ]] || die "--mem-gb $MEM_GB is more than this machine has ($HOST_MEM_GB GB)"

        build_mounts
        echo "Image           : ${IMAGE}  (code from ${MODE}: ${CODE_HOME})"
        echo "Container limits: ${CORES} cores, ${MEM_GB} GB RAM (hard)"
        echo "Mounted folders : $(printf '%s ' "${VOLUME_ARGS[@]}" | sed 's/-v //g; s/:[^ ]*//g')"
        docker_run "$CORES" "$MEM_GB" python "${CODE_HOME}/workflow/scripts/run_batches.py" \
            --cores "$CORES" --mem-gb "$MEM_GB" "${RUN_ARGS[@]}"
        ;;

    *)
        die "unknown command '$COMMAND'. Use: samplesheet, run, shell, build (or --help)"
        ;;
esac
