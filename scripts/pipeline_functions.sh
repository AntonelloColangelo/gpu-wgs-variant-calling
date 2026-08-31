#!/usr/bin/env bash
# Technical functions used by run_parabricks_hg002.sh.
# The main file can then be read on its own: Docker, checkpoints, logging and
# safe export all live here.

set -Eeuo pipefail

say() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$*"
}

die() {
    printf '\nERROR: %s\n' "$*" >&2
    exit 1
}


# -----------------------------------------------------------------------------
# Checkpoint provenance
#
# Checkpoints are indexed on PREFIX and WORK_VOLUME. If those names do not depend
# on the inputs, changing the reference and re-launching the pipeline finds the
# old BAM again, which passes validation because it is structurally sound:
# alignment gets skipped and the pipeline reports PASS while returning exactly
# the results it was supposed to replace.
#
# The key computed here ties prefix and volumes to everything that determines the
# result. Change one ingredient and you get new names: the old checkpoints are no
# longer found and at the same time are not destroyed, so the previous run stays
# available as a baseline.
# -----------------------------------------------------------------------------

require_runtime() {
    grep -qi microsoft /proc/version ||
        die "run this script under Ubuntu WSL2, not under Git Bash."
    command -v docker >/dev/null || die "Docker is not available in Ubuntu WSL2."
    docker info >/dev/null 2>&1 || die "Docker Desktop is not responding."
    docker image inspect "$PB_IMAGE" >/dev/null 2>&1 ||
        die "Parabricks image not found: $PB_IMAGE"
}

# One line per ingredient, name and value separated by a tab. The same function
# feeds both the key computation and the plain-text manifest: a hash on its own
# is not verifiable, you need to be able to read what went into it.
#
# The .fai is the central choice: it weighs a few kB and changes if even a single
# contig of the reference changes, so it identifies the reference without having
# to read its 3 Gb.
#
# WARNING: every literal below is hashed. The Italian fallback 'nessuno' stays as
# it is on purpose — translating it would change the key of every run without
# intervals, which is how the published run 53007e55 and its 92 GB volume are
# addressed. The manifest has to keep saying what was actually hashed.
run_key_ingredients() {
    printf 'parabricks_image_id\t%s\n' \
        "$(docker image inspect --format '{{.Id}}' "$PB_IMAGE")"
    printf 'reference\t%s\n' "$(basename "$REFERENCE_HOST")"
    printf 'reference_fai_sha256\t%s\n' \
        "$(sha256sum "${REFERENCE_HOST}.fai" | cut -d' ' -f1)"
    printf 'known_sites\t%s %s\n' \
        "$(basename "$KNOWN_SITES_HOST")" "$(stat -c '%s' "$KNOWN_SITES_HOST")"
    printf 'fastq_r1\t%s %s\n' \
        "$(basename "$FASTQ_R1_HOST")" "$(stat -c '%s' "$FASTQ_R1_HOST")"
    printf 'fastq_r2\t%s %s\n' \
        "$(basename "$FASTQ_R2_HOST")" "$(stat -c '%s' "$FASTQ_R2_HOST")"
    printf 'read_group\t%s\n' "$READ_GROUP"
    printf 'intervals\t%s\n' "${INTERVAL_ARGS[*]:-nessuno}"
}

compute_run_key() {
    local f
    for f in "$REFERENCE_HOST" "${REFERENCE_HOST}.fai" \
             "$KNOWN_SITES_HOST" "$FASTQ_R1_HOST" "$FASTQ_R2_HOST"; do
        [[ -s "$f" ]] || die "a file needed to identify the run is missing: $f
   If it is the reference, download it with: bash scripts/fetch_reference_noalt.sh"
    done
    run_key_ingredients | sha256sum | cut -c1-8
}

json_escape() {
    local s="$1"
    s="${s//\\/\\\\}"
    s="${s//\"/\\\"}"
    printf '%s' "$s"
}

write_run_manifest() {
    local path="$1"
    local first=1 key value
    {
        printf '{\n'
        printf '  "run_id": "%s",\n'      "$(json_escape "$RUN_ID")"
        printf '  "mode": "%s",\n'        "$(json_escape "$MODE")"
        printf '  "run_key": "%s",\n'     "$(json_escape "$RUN_KEY")"
        printf '  "prefix": "%s",\n'      "$(json_escape "$PREFIX")"
        printf '  "work_volume": "%s",\n' "$(json_escape "$WORK_VOLUME")"
        printf '  "ingredients": {\n'
        while IFS=$'\t' read -r key value; do
            [[ -n "$key" ]] || continue
            (( first )) || printf ',\n'
            first=0
            printf '    "%s": "%s"' "$(json_escape "$key")" "$(json_escape "$value")"
        done < <(run_key_ingredients)
        printf '\n  }\n}\n'
    } > "${path}.tmp"
    mv -f "${path}.tmp" "$path"
}

write_state() {
    local key="$1"
    local label="$2"
    local status="$3"
    local started="${4:-}"
    printf 'STEP_KEY=%s\nSTEP_LABEL=%s\nSTEP_STATUS=%s\nSTEP_STARTED=%s\n' \
        "$key" "$label" "$status" "$started" > "${STATE_FILE}.tmp"
    mv -f "${STATE_FILE}.tmp" "${STATE_FILE}"
}

record_timing() {
    local key="$1"
    local started="$2"
    local ended="$3"
    local seconds="$4"
    local status="$5"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$key" "$started" "$ended" "$seconds" "$status" >> "$TIMING_FILE"
}

skip_step() {
    local key="$1"
    local label="$2"
    local now
    now="$(date --iso-8601=seconds)"
    printf '[OK] %s: checkpoint already valid, moving on.\n' "$label"
    record_timing "$key" "$now" "$now" "0" "SKIPPED_VALID"
}

run_step() {
    local key="$1"
    local label="$2"
    shift 2

    local started ended start_epoch end_epoch rc log_file
    started="$(date --iso-8601=seconds)"
    start_epoch="$(date +%s)"
    log_file="${RUN_DIR}/${key}.log"

    say "$label"
    write_state "$key" "$label" "RUNNING" "$started"

    set +e
    "$@" 2>&1 | tee "$log_file"
    rc="${PIPESTATUS[0]}"
    set -e

    ended="$(date --iso-8601=seconds)"
    end_epoch="$(date +%s)"

    if [[ "$rc" -eq 0 ]]; then
        record_timing "$key" "$started" "$ended" "$((end_epoch - start_epoch))" "PASS"
        write_state "$key" "$label" "PASS" "$started"
        printf '[OK] %s complete.\n' "$label"
        return 0
    fi

    record_timing "$key" "$started" "$ended" "$((end_epoch - start_epoch))" "FAILED"
    write_state "$key" "$label" "FAILED" "$started"
    printf '[ERROR] %s failed (code %s). Log: %s\n' "$label" "$rc" "$log_file" >&2
    return "$rc"
}

remove_stopped_container() {
    local name="$1"
    local running
    running="$(docker inspect -f '{{.State.Running}}' "$name" 2>/dev/null || true)"
    if [[ "$running" == "true" ]]; then
        die "container ${name} is already running. I will not start a second process on the same files."
    fi
    if [[ "$running" == "false" ]]; then
        docker rm "$name" >/dev/null
    fi
}

run_parabricks() {
    local container_name="$1"
    shift
    remove_stopped_container "$container_name"
    docker run --rm \
        --name "$container_name" \
        --gpus all \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${WORK_VOLUME},target=/work" \
        --mount "type=volume,source=${TMP_VOLUME},target=/pbtmp" \
        "$PB_IMAGE" "$@"
}

# postprocess.sh and generate_report.py receive the reference from the
# environment instead of writing it inside themselves: changing it in one place
# is enough for the whole pipeline.
run_tools() {
    local container_name="$1"
    shift
    remove_stopped_container "$container_name"
    docker run --rm \
        --name "$container_name" \
        --env "REFERENCE=${REFERENCE}" \
        --env "REFERENCE_DICT=${REFERENCE_DICT}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${WORK_VOLUME},target=/work" \
        --mount "type=volume,source=${TMP_VOLUME},target=/pbtmp" \
        "$TOOLS_IMAGE" "$@"
}

tools_quiet() {
    docker run --rm \
        --env "REFERENCE=${REFERENCE}" \
        --env "REFERENCE_DICT=${REFERENCE_DICT}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${WORK_VOLUME},target=/work" \
        --mount "type=volume,source=${TMP_VOLUME},target=/pbtmp" \
        "$TOOLS_IMAGE" "$@" >/dev/null 2>&1
}

work_remove() {
    if [[ "$#" -gt 0 ]]; then
        docker run --rm \
            --mount "type=volume,source=${WORK_VOLUME},target=/work" \
            alpine:3.20 rm -f -- "$@" >/dev/null
    fi
}

work_file_exists() {
    local path="$1"
    docker run --rm \
        --mount "type=volume,source=${WORK_VOLUME},target=/work,readonly" \
        alpine:3.20 test -s "$path" >/dev/null 2>&1
}

valid_fastq_json() {
    local path="/work/${PREFIX}.fastq_validation.json"
    tools_quiet python3 -c '
import json, sys
d = json.load(open(sys.argv[1], encoding="utf-8"))
expected = int(sys.argv[2])
ok = (
    d.get("status") == "PASS"
    and d.get("full_gzip_integrity_checked") is True
    and d.get("paired_ids_match") is True
    and d.get("fastq_structure_valid") is True
    and int(d.get("total_pairs", -1)) == expected
)
raise SystemExit(0 if ok else 1)
' "$path" "$EXPECTED_PAIRS"
}

valid_bam() {
    local bam="/work/${PREFIX}.bam"
    local bai="${bam}.bai"
    tools_quiet bash /project/scripts/postprocess.sh check-bam "$bam" "$bai" "$SAMPLE" "$READ_GROUP_PU"
}

valid_partial_bam() {
    local bam="/work/${PREFIX}.partial.bam"
    local bai="${bam}.bai"
    tools_quiet bash /project/scripts/postprocess.sh check-bam "$bam" "$bai" "$SAMPLE" "$READ_GROUP_PU"
}

valid_recal() {
    tools_quiet bash /project/scripts/postprocess.sh check-recal "/work/${PREFIX}.recal.txt"
}

valid_vcf() {
    local path="$1"
    tools_quiet bash /project/scripts/postprocess.sh check-vcf "$path" "${path}.tbi" "$SAMPLE"
}

valid_hardfilter() {
    tools_quiet bash /project/scripts/postprocess.sh check-hardfilter "$PREFIX" "$SAMPLE"
}

valid_annotation() {
    tools_quiet bash /project/scripts/postprocess.sh check-annotation "$PREFIX" "$SAMPLE"
}

valid_qc() {
    tools_quiet bash /project/scripts/postprocess.sh check-qc "$PREFIX"
}

copy_validation_checkpoint() {
    local candidate relative
    for candidate in "${VALIDATION_CANDIDATES[@]}"; do
        [[ -s "$candidate" ]] || continue
        relative="${candidate#"$PROJECT_DIR"/}"
        docker run --rm \
            --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
            --mount "type=volume,source=${WORK_VOLUME},target=/work" \
            alpine:3.20 cp "/project/${relative}" "/work/${PREFIX}.fastq_validation.json"
        if valid_fastq_json; then
            printf 'Reusing the full check that already passed: %s\n' "$relative"
            return 0
        fi
        work_remove "/work/${PREFIX}.fastq_validation.json"
    done
    return 1
}

ensure_fastq_integrity() {
    if valid_fastq_json; then
        skip_step "1_fastq_integrity" "1/9 - FASTQ integrity and pairing"
        return
    fi

    if [[ "$(stat -c '%s' "$FASTQ_R1_HOST")" == "$EXPECTED_R1_BYTES" ]] &&
       [[ "$(stat -c '%s' "$FASTQ_R2_HOST")" == "$EXPECTED_R2_BYTES" ]] &&
       copy_validation_checkpoint; then
        skip_step "1_fastq_integrity" "1/9 - FASTQ integrity and pairing"
        return
    fi

    work_remove "/work/${PREFIX}.fastq_validation.json"
    run_step "1_fastq_integrity" "1/9 - FASTQ integrity and pairing" \
        run_tools "hg002_1_fastq_integrity" \
        bash /project/scripts/postprocess.sh integrity \
        "$PREFIX" "$FASTQ_R1_CONTAINER" "$FASTQ_R2_CONTAINER"

    valid_fastq_json || die "the FASTQ check did not produce a valid checkpoint."
}

# The check iteration 1 lacked, and which would have cost zero seconds.
#
# A reference with _alt contigs but without the .alt file beside it makes BWA
# treat ALT and primary locus as ordinary multi-mapping: MAPQ 0, and
# HaplotypeCaller discards those reads below MAPQ 20. In the first iteration the
# defect only surfaced after the benchmark, ten hours later, and it cost 98.51 %
# of the true SNPs in the MHC of chr6.
#
# To reproduce the defective baseline deliberately:
#   HG002_ALLOW_ALT_WITHOUT_ALT_FILE=1 HG002_REF_NAME=Homo_sapiens_assembly38.fasta ...
check_reference_alt_contigs() {
    local total n_alt
    total="$(wc -l < "${REFERENCE_HOST}.fai")"
    n_alt="$(cut -f1 "${REFERENCE_HOST}.fai" | grep -c '_alt$' || true)"

    printf 'Reference: %s\n' "$(basename "$REFERENCE_HOST")"
    printf '  total contigs: %s | _alt contigs: %s\n' "$total" "$n_alt"

    if (( n_alt == 0 )); then
        echo "  no _alt contigs: alt-aware alignment is not needed."
        return 0
    fi
    if [[ -s "${REFERENCE_HOST}.alt" ]]; then
        echo "  .alt file present: BWA can work alt-aware."
        return 0
    fi
    if [[ "${HG002_ALLOW_ALT_WITHOUT_ALT_FILE:-0}" == "1" ]]; then
        printf '  WARNING: %s _alt contigs without an .alt file. Continuing because it was\n' "$n_alt"
        printf '  requested explicitly: ambiguous reads will get MAPQ 0 and the variants\n'
        printf '  underneath will not be called.\n'
        return 0
    fi
    die "the reference contains ${n_alt} _alt contigs but $(basename "$REFERENCE_HOST").alt is missing.
   Without that file BWA assigns MAPQ 0 to reads mapping both to the primary
   locus and to an alternate contig, and HaplotypeCaller discards them: in the
   first iteration that cost 98.51 % of the SNPs in the MHC of chr6.
   Fix: use a no-ALT reference (bash scripts/fetch_reference_noalt.sh) or obtain
   the .alt file.
   To reproduce the defective baseline on purpose, re-run with
   HG002_ALLOW_ALT_WITHOUT_ALT_FILE=1."
}

preflight() {
    say "Preflight checks"

    require_runtime
    command -v nvidia-smi >/dev/null || die "nvidia-smi is not available in Ubuntu WSL2."

    # The list is derived from REFERENCE_HOST: changing reference does not mean
    # updating ten paths by hand.
    local required=(
        "$FASTQ_R1_HOST"
        "$FASTQ_R2_HOST"
        "$REFERENCE_HOST"
        "${REFERENCE_HOST}.fai"
        "$REFERENCE_DICT_HOST"
        "${REFERENCE_HOST}.amb"
        "${REFERENCE_HOST}.ann"
        "${REFERENCE_HOST}.bwt"
        "${REFERENCE_HOST}.pac"
        "${REFERENCE_HOST}.sa"
        "$KNOWN_SITES_HOST"
        "${KNOWN_SITES_HOST}.tbi"
        "$PROJECT_DIR/snpEff_data/hg38/snpEffectPredictor.bin"
        "$PROJECT_DIR/scripts/postprocess.sh"
        "$PROJECT_DIR/scripts/generate_report.py"
    )
    local file
    for file in "${required[@]}"; do
        [[ -s "$file" ]] || die "required file missing or empty: $file"
    done

    check_reference_alt_contigs

    docker image inspect "$TOOLS_IMAGE" >/dev/null ||
        die "bioinformatics image not found: $TOOLS_IMAGE"
    docker image inspect alpine:3.20 >/dev/null ||
        die "utility image alpine:3.20 not found."

    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader |
        sed 's/^/GPU: /'

    docker volume create "$WORK_VOLUME" >/dev/null
    docker volume create "$TMP_VOLUME" >/dev/null

    local docker_free_kb host_free_kb docker_min_kb host_min_kb
    docker_free_kb="$(
        docker run --rm \
            --mount "type=volume,source=${WORK_VOLUME},target=/work" \
            alpine:3.20 sh -c "df -Pk /work | awk 'END {print \$4}'"
    )"
    host_free_kb="$(df -Pk "$PROJECT_DIR" | awk 'END {print $4}')"

    if [[ "$MODE" == "full" ]]; then
        docker_min_kb=$((300 * 1024 * 1024))
        host_min_kb=$((150 * 1024 * 1024))
    else
        docker_min_kb=$((5 * 1024 * 1024))
        host_min_kb=$((5 * 1024 * 1024))
    fi

    (( docker_free_kb >= docker_min_kb )) ||
        die "not enough space on Docker's Linux disk."
    (( host_free_kb >= host_min_kb )) ||
        die "not enough space on the Windows disk to export the results."

    printf 'Free space, Docker: %.1f GiB\n' "$(awk -v k="$docker_free_kb" 'BEGIN{print k/1024/1024}')"
    printf 'Free space, Windows: %.1f GiB\n' "$(awk -v k="$host_free_kb" 'BEGIN{print k/1024/1024}')"
    printf 'Mode: %s | Sample: %s | Results volume: %s\n' \
        "$MODE" "$SAMPLE" "$WORK_VOLUME"
    printf 'Provenance key: %s | Prefix: %s\n' "$RUN_KEY" "$PREFIX"
}

finalize_bam() {
    run_tools "hg002_2_bam_check" \
        bash /project/scripts/postprocess.sh finalize-bam \
        "/work/${PREFIX}.partial.bam" \
        "/work/${PREFIX}.bam" \
        "$SAMPLE" "$READ_GROUP_PU" \
        "/work/${PREFIX}.duplicate_metrics.partial.txt" \
        "/work/${PREFIX}.duplicate_metrics.txt"
}

finalize_recal() {
    run_tools "hg002_3_bqsr_check" \
        bash /project/scripts/postprocess.sh finalize-recal \
        "/work/${PREFIX}.recal.partial.txt" \
        "/work/${PREFIX}.recal.txt"
}

finalize_vcf() {
    local partial="$1"
    local final="$2"
    local name="$3"
    run_tools "$name" \
        bash /project/scripts/postprocess.sh finalize-vcf \
        "$partial" "$final" "$SAMPLE"
}

export_results() {
    local stage="${OUTPUT_DIR}/.incoming_${PREFIX}_${RUN_ID}"
    mkdir -p "$stage"

    local files=(
        "${PREFIX}.fastq_validation.json"
        "${PREFIX}.bam"
        "${PREFIX}.bam.bai"
        "${PREFIX}.duplicate_metrics.txt"
        "${PREFIX}.recal.txt"
        "${PREFIX}.g.vcf.gz"
        "${PREFIX}.g.vcf.gz.tbi"
        "${PREFIX}.vcf.gz"
        "${PREFIX}.vcf.gz.tbi"
        "${PREFIX}.hardfiltered.vcf.gz"
        "${PREFIX}.hardfiltered.vcf.gz.tbi"
        "${PREFIX}_annotated.vcf.gz"
        "${PREFIX}_annotated.vcf.gz.tbi"
        "${PREFIX}_high_impact_variants.tsv"
        "${PREFIX}_snpEff_summary.html"
        "${PREFIX}_snpEff_summary.genes.txt"
        "${PREFIX}.snpeff.log"
        "${PREFIX}.annotation_qc.tsv"
        "${PREFIX}.flagstat.txt"
        "${PREFIX}.samtools_stats.txt"
        "${PREFIX}.idxstats.tsv"
        "${PREFIX}.bam_header.sam"
        "${PREFIX}.coverage.tsv"
        "${PREFIX}.wgs_metrics.txt"
        "${PREFIX}.depth_thresholds.tsv"
        "${PREFIX}.vcf_contigs.txt"
    )

    local source_paths=()
    local name
    for name in "${files[@]}"; do
        source_paths+=("/work/${name}")
    done

    docker run --rm \
        --mount "type=volume,source=${WORK_VOLUME},target=/work,readonly" \
        --mount "type=bind,source=${stage},target=/export" \
        alpine:3.20 cp -a "${source_paths[@]}" /export/

    cp "$TIMING_FILE" "${stage}/${PREFIX}.step_times.tsv"

    docker run --rm \
        --mount "type=bind,source=${stage},target=/export,readonly" \
        "$TOOLS_IMAGE" bash -c '
set -euo pipefail
p="$1"
samtools quickcheck -v "/export/${p}.bam"
test -s "/export/${p}.bam.bai"
samtools idxstats "/export/${p}.bam" >/dev/null
for v in \
    "/export/${p}.g.vcf.gz" \
    "/export/${p}.vcf.gz" \
    "/export/${p}.hardfiltered.vcf.gz" \
    "/export/${p}_annotated.vcf.gz"
do
    gzip -t "$v"
    test -s "${v}.tbi"
    tabix -l "$v" >/dev/null
done
' _ "$PREFIX"

    for name in "${files[@]}" "${PREFIX}.step_times.tsv"; do
        mv -f "${stage}/${name}" "${OUTPUT_DIR}/${name}"
    done
    rmdir "$stage"
}

generate_report() {
    docker run --rm \
        --name "hg002_9_report" \
        --entrypoint python3 \
        --env "REFERENCE=${REFERENCE}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=bind,source=${OUTPUT_DIR},target=/project/output" \
        --mount "type=bind,source=${REPORT_DIR},target=/project/reports" \
        "$TOOLS_IMAGE" \
        /project/scripts/generate_report.py report \
        --root /project \
        --prefix "$PREFIX" \
        --sample "$SAMPLE" \
        --mode "$MODE"
}

pipeline_failed() {
    local rc="$1"
    local line="$2"
    printf '\nPipeline aborted at line %s (code %s).\n' "$line" "$rc" >&2
    printf 'The results already validated stay in volume %s and will be reused.\n' \
        "$WORK_VOLUME" >&2
}
