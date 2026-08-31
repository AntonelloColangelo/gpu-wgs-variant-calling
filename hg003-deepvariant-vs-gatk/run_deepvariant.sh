#!/usr/bin/env bash
# =============================================================================
# DeepVariant on the BAM the pipeline has already produced.
#
# This is an experiment, not a stage: it answers "is HaplotypeCaller's F1 the
# ceiling of this data, or only the ceiling of this caller?" and it answers it
# by measurement rather than by citing someone else's benchmark.
#
# It reuses the existing alignment. That is the whole reason this is cheap: the
# 82 GB BAM in the work volume already exists, duplicates marked, and with the
# BQSR table applied only in memory during calling -- so on disk it carries no
# recalibrated qualities, which is exactly the input DeepVariant wants. Nothing
# is realigned, nothing the pipeline produced is touched or overwritten.
#
#   Usage:  bash scripts/run_deepvariant.sh
#           PREFIX=... bash scripts/run_deepvariant.sh
#
# Then compare the two callers on the same truth set:
#
#   QUERY_VCF=output/<prefix>.deepvariant.vcf.gz BENCH_LABEL=<prefix>_deepvariant \
#       bash benchmark_giab.sh
# =============================================================================

set -Eeuo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${PROJECT_DIR}/config/sample.env}"
[[ -s "$CONFIG_FILE" ]] || { echo "Configuration not found: $CONFIG_FILE" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG_FILE"
: "${SAMPLE_ID:?SAMPLE_ID is required in the configuration}"

REF_NAME="${PARABRICKS_REF_NAME:-${REFERENCE_NAME:-GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}}"
REFERENCE_HOST="${PROJECT_DIR}/ref/${REF_NAME}"
REFERENCE="/project/ref/${REF_NAME}"
PB_IMAGE="${PB_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1}"
SAMPLE_SLUG="${SAMPLE_ID,,}"

msg() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31m!! %s\033[0m\n' "$*" >&2; exit 1; }

# --- 1. Which run to call on -------------------------------------------------
# The prefix carries the provenance key of the alignment, so it also identifies
# the volume the BAM lives in. Inferred when unambiguous, never guessed.
if [[ -z "${PREFIX:-}" ]]; then
    mapfile -t candidates < <(
        find "${PROJECT_DIR}/output" -maxdepth 1 -name '*.bam' -printf '%f\n' 2>/dev/null |
            sed 's/\.bam$//' | grep -v 'smoke' | sort
    )
    case "${#candidates[@]}" in
        0) die "No BAM in output/. Run the pipeline first." ;;
        1) PREFIX="${candidates[0]}" ;;
        *) printf 'More than one BAM available: pick one with PREFIX=...\n' >&2
           printf '  %s\n' "${candidates[@]}" >&2
           exit 1 ;;
    esac
fi

RUN_KEY="${PREFIX##*_}"
WORK_VOLUME="${WORK_VOLUME:-${SAMPLE_SLUG}_work_${RUN_KEY}}"
TMP_VOLUME="${TMP_VOLUME:-${SAMPLE_SLUG}_dv_tmp_${RUN_KEY}}"
OUTPUT_DIR="${PROJECT_DIR}/output"
LOG_DIR="${PROJECT_DIR}/logs/deepvariant"
LOG_FILE="${LOG_DIR}/${PREFIX}.deepvariant.log"

FINAL_VCF="/work/${PREFIX}.deepvariant.vcf.gz"
PARTIAL_VCF="/work/${PREFIX}.deepvariant.partial.vcf.gz"

mkdir -p "$LOG_DIR"

# Questo lavoro dura ore e gira dentro un container: senza uno stato scritto e
# senza un comando da seguire, dal terminale è indistinguibile da un processo
# bloccato. Scrive lo stesso file che leggono Get-PipelineStatus.ps1 e
# Watch-Parabricks.ps1, così il riscontro visivo non dipende da chi lo lancia.
STATE_FILE="${PROJECT_DIR}/logs/current_step.env"
write_state() {
    printf 'STEP_KEY=%s\nSTEP_LABEL=%s\nSTEP_STATUS=%s\nSTEP_STARTED=%s\n' \
        "$1" "$2" "$3" "$(date --iso-8601=seconds)" > "${STATE_FILE}.tmp"
    mv -f "${STATE_FILE}.tmp" "$STATE_FILE"
}
on_exit() {
    local rc="$?"
    [[ "$rc" -eq 0 ]] || write_state "deepvariant" "DeepVariant su ${PREFIX}" "FAILED"
    return "$rc"
}
trap on_exit EXIT

# Lo stesso lock della pipeline: entrambe vogliono la GPU, e due job insieme non
# ci stanno. Senza, si sovrascriverebbero anche lo stato qui sopra.
exec 9>"${PROJECT_DIR}/logs/pipeline.lock"
flock -n 9 || die "un altro lavoro di questo progetto è già in corso."

# --- 2. Checks ---------------------------------------------------------------
msg "1/4 - Checks"
grep -qi microsoft /proc/version || die "run this under Ubuntu WSL2, not Git Bash."
command -v docker >/dev/null || die "Docker is not available."
docker info >/dev/null 2>&1 || die "Docker Desktop is not responding."
docker image inspect "$PB_IMAGE" >/dev/null 2>&1 || die "Parabricks image missing: $PB_IMAGE"
[[ -s "$REFERENCE_HOST" ]] || die "Reference not found: $REFERENCE_HOST"
docker volume inspect "$WORK_VOLUME" >/dev/null 2>&1 ||
    die "Work volume not found: $WORK_VOLUME"

docker run --rm --mount "type=volume,source=${WORK_VOLUME},target=/work,readonly" \
    alpine:3.20 test -s "/work/${PREFIX}.bam" ||
    die "BAM not found in volume ${WORK_VOLUME}: ${PREFIX}.bam"

echo "Prefix:    ${PREFIX}"
echo "Volume:    ${WORK_VOLUME}"
echo "Reference: ${REF_NAME}"
cat <<EOF

Per seguirlo dal terminale VSCode, in PowerShell dalla cartella del progetto:

  Get-Content -Wait -Tail 20 "logs\\deepvariant\\${PREFIX}.deepvariant.log"

oppure, per lo stato sintetico:  .\\Get-PipelineStatus.ps1 -Watch
EOF

if docker run --rm --mount "type=volume,source=${WORK_VOLUME},target=/work,readonly" \
        alpine:3.20 test -s "/work/${PREFIX}.deepvariant.vcf.gz"; then
    echo "DeepVariant VCF already present in the volume: nothing to call."
    SKIP_CALLING=1
else
    SKIP_CALLING=0
fi

# --- 3. Calling --------------------------------------------------------------
# Deliberately NOT using --enable-small-model. It routes easy sites through a
# lighter network and is the one real speed lever here, but this run exists to
# measure accuracy: the honest comparison against HaplotypeCaller is the default
# model. Turn it on later, once you know what the ceiling costs in time.
#
# No BQSR table is passed, and that is not an omission. DeepVariant learns the
# behaviour of base qualities from the pileup itself; feeding it recalibrated
# qualities is documented as making its results worse, not better.
if [[ "$SKIP_CALLING" -eq 0 ]]; then
    msg "2/4 - pbrun deepvariant (hours, on one GPU)"
    write_state "deepvariant" "DeepVariant su ${PREFIX}" "RUNNING"
    docker rm -f "${SAMPLE_SLUG}_deepvariant" >/dev/null 2>&1 || true
    docker run --rm \
        --name "${SAMPLE_SLUG}_deepvariant" \
        --gpus all \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${WORK_VOLUME},target=/work" \
        --mount "type=volume,source=${TMP_VOLUME},target=/pbtmp" \
        "$PB_IMAGE" \
        pbrun deepvariant \
            --ref "$REFERENCE" \
            --in-bam "/work/${PREFIX}.bam" \
            --out-variants "$PARTIAL_VCF" \
            --tmp-dir /pbtmp \
            --num-gpus 1 \
        2>&1 | tee "$LOG_FILE"

    msg "3/4 - Check and finalisation"
    write_state "deepvariant_check" "Verifica del VCF DeepVariant" "RUNNING"
    # Promoted to its final name only after the file has been proved readable
    # end to end: a truncated VCF that gzip accepts is exactly the failure the
    # rest of this project is built to prevent.
    docker run --rm \
        --mount "type=volume,source=${WORK_VOLUME},target=/work" \
        --env "SAMPLE=${SAMPLE_ID}" \
        "${TOOLS_IMAGE:-bioinfo-codeserver:latest}" bash -euc '
p="$1"
cd /work
gzip -t "${p}.deepvariant.partial.vcf.gz"
tabix -f -p vcf "${p}.deepvariant.partial.vcf.gz"
tabix -l "${p}.deepvariant.partial.vcf.gz" >/dev/null
tabix -H "${p}.deepvariant.partial.vcf.gz" | awk -F "\t" -v s="$SAMPLE" "
    /^#CHROM/ { found=1; if (NF < 10 || \$10 != s) exit 2 }
    END { if (!found) exit 3 }"
n=$(zcat "${p}.deepvariant.partial.vcf.gz" | grep -vc "^#")
[ "$n" -gt 0 ] || { echo "empty VCF" >&2; exit 2; }
printf "records: %s\n" "$n"
mv -f "${p}.deepvariant.partial.vcf.gz"     "${p}.deepvariant.vcf.gz"
mv -f "${p}.deepvariant.partial.vcf.gz.tbi" "${p}.deepvariant.vcf.gz.tbi"
' _ "$PREFIX" || die "the DeepVariant VCF did not pass its checks."
else
    msg "2-3/4 - Calling skipped, result already validated"
fi

# --- 4. Export ---------------------------------------------------------------
msg "4/4 - Export to Windows"
STAGE_DIR="${OUTPUT_DIR}/.incoming_dv_${PREFIX}"
mkdir -p "$STAGE_DIR"
docker run --rm \
    --mount "type=volume,source=${WORK_VOLUME},target=/work,readonly" \
    --mount "type=bind,source=${STAGE_DIR},target=/export" \
    alpine:3.20 cp -a \
        "/work/${PREFIX}.deepvariant.vcf.gz" \
        "/work/${PREFIX}.deepvariant.vcf.gz.tbi" /export/
mv -f "${STAGE_DIR}/${PREFIX}.deepvariant.vcf.gz"     "${OUTPUT_DIR}/"
mv -f "${STAGE_DIR}/${PREFIX}.deepvariant.vcf.gz.tbi" "${OUTPUT_DIR}/"
rmdir "$STAGE_DIR"

docker volume rm "$TMP_VOLUME" >/dev/null 2>&1 || true
write_state "deepvariant_complete" "DeepVariant completato su ${PREFIX}" "PASS"

msg "Done"
cat <<EOF
VCF: output/${PREFIX}.deepvariant.vcf.gz
Log: ${LOG_FILE#"${PROJECT_DIR}/"}

Now measure it against the same truth set HaplotypeCaller was measured on:

  QUERY_VCF=output/${PREFIX}.deepvariant.vcf.gz \\
  BENCH_LABEL=${PREFIX}_deepvariant \\
      bash benchmark_giab.sh
EOF
