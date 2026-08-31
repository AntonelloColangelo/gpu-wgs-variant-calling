#!/usr/bin/env bash
# =============================================================================
# What the GPU is worth on this machine: same hardware, same data, same
# reference, one variable only (what does the computing).
#
#   Test A - alignment:      pbrun fq2bam            vs  bwa mem + sort + MarkDuplicates
#   Test B - variant calling: pbrun haplotypecaller  vs  gatk HaplotypeCaller
#
# Test A runs on a subset of read pairs: the ratio between the two timings does
# not depend on how many reads there are, while a full CPU run would cost more
# than a day of machine time.
#
# Test B runs on the full BAM the pipeline already produced, so at real coverage
# (~35x) and not at the thin coverage a subset would give. Neither caller
# applies BQSR: what is compared is the caller, not the recalibration.
#
# Usage:
#   bash scripts/benchmark_cpu_vs_gpu.sh
#   PAIRS=5000000 REGION=chr20:1-5000000 bash scripts/benchmark_cpu_vs_gpu.sh
# =============================================================================

set -Eeuo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"

# The defaults are the ones behind the published numbers: change them and the
# results stop being comparable with the table in the README.
PAIRS="${PAIRS:-10000000}"
REGION="${REGION:-chr20:1-20000000}"
THREADS="${THREADS:-$(nproc)}"

REF_NAME="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
REFERENCE="/project/ref/${REF_NAME}"

PB_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
TOOLS_IMAGE="bioinfo-codeserver:latest"

# The full-coverage BAM of iteration 2, already validated by the pipeline.
FULL_VOLUME="hg002_work_53007e55"
FULL_BAM="/full/HG002_NovaSeq_40x_53007e55.bam"

BENCH_VOLUME="${BENCH_VOLUME:-hg002_bench}"
BENCH_TMP_VOLUME="${BENCH_TMP_VOLUME:-hg002_bench_tmp}"

READ_GROUP='@RG\tID:HV3C3DSXX.2\tPL:ILLUMINA\tPM:NovaSeq6000\tLB:HG002_PCR_FREE\tPU:HV3C3DSXX.2.AGCGATAG+AGGCGAAG\tSM:HG002'

FASTQ_R1="/project/HG002.novaseq.pcr-free.40x.R1.fastq.gz"
FASTQ_R2="/project/HG002.novaseq.pcr-free.40x.R2.fastq.gz"

RESULT_DIR="${PROJECT_DIR}/reports/cpu_vs_gpu"
LOG_DIR="${PROJECT_DIR}/logs/bench"
RESULT_TSV="${RESULT_DIR}/benchmark.tsv"
RESULT_JSON="${RESULT_DIR}/benchmark.json"

mkdir -p "$RESULT_DIR" "$LOG_DIR"


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# Timings end up in a single file, one line per phase, so that the figure and
# the README read the same source.
record() {
    local phase="$1" seconds="$2" note="$3"
    printf '%s\t%s\t%s\n' "$phase" "$seconds" "$note" >> "$RESULT_TSV"
    log "  -> ${phase}: ${seconds}s (${note})"
}

# Times a command, keeps its log and returns the seconds in TIMED.
TIMED=0
timed() {
    local logfile="$1"
    shift
    local start end
    start="$(date +%s)"
    if ! "$@" > "$logfile" 2>&1; then
        printf 'command failed, last lines of %s:\n' "$logfile" >&2
        tail -25 "$logfile" >&2
        die "phase aborted"
    fi
    end="$(date +%s)"
    TIMED="$((end - start))"
}

on_bench() {
    docker run --rm \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${BENCH_VOLUME},target=/bench" \
        "$TOOLS_IMAGE" "$@"
}


# -----------------------------------------------------------------------------
# Preflight
# -----------------------------------------------------------------------------

log "CPU vs GPU benchmark - ${PAIRS} pairs, region ${REGION}, ${THREADS} threads"

docker image inspect "$PB_IMAGE" >/dev/null 2>&1 || die "image missing: $PB_IMAGE"
docker image inspect "$TOOLS_IMAGE" >/dev/null 2>&1 || die "image missing: $TOOLS_IMAGE"
docker volume inspect "$FULL_VOLUME" >/dev/null 2>&1 ||
    die "volume missing: $FULL_VOLUME (test B needs the full BAM)"
[[ -s "${PROJECT_DIR}/ref/${REF_NAME}" ]] || die "reference missing: ref/${REF_NAME}"
[[ -s "${PROJECT_DIR}/ref/${REF_NAME}.bwt" ]] || die "BWA index missing next to the reference"

docker volume create "$BENCH_VOLUME" >/dev/null
docker volume create "$BENCH_TMP_VOLUME" >/dev/null

printf 'phase\tseconds\tnote\n' > "$RESULT_TSV"


# -----------------------------------------------------------------------------
# Preparation: deterministic subset of the FASTQ files
#
# The first PAIRS pairs of both files. head closes the pipe and kills zcat with
# SIGPIPE: that is expected, so pipefail stays off here.
# -----------------------------------------------------------------------------

SUB_R1="/bench/sub_R1.fastq.gz"
SUB_R2="/bench/sub_R2.fastq.gz"
LINES="$((PAIRS * 4))"

if on_bench test -s "$SUB_R1" && on_bench test -s "$SUB_R2"; then
    log "Subset already present, reusing it"
else
    log "Extracting ${PAIRS} pairs from the full FASTQ files"
    timed "${LOG_DIR}/0_subset.log" \
        on_bench bash -c "set +o pipefail
            zcat '${FASTQ_R1}' | head -n ${LINES} | gzip -1 > '${SUB_R1}'
            zcat '${FASTQ_R2}' | head -n ${LINES} | gzip -1 > '${SUB_R2}'"
    record "prep_subset" "$TIMED" "extraction of ${PAIRS} pairs"
fi


# -----------------------------------------------------------------------------
# TEST A1 - alignment on GPU: pbrun fq2bam
#
# Same flags as the production pipeline, including the limits imposed by RAM:
# that is the configuration whose cost we want to know.
# -----------------------------------------------------------------------------

log "Test A1 - pbrun fq2bam (GPU)"
timed "${LOG_DIR}/A1_gpu_fq2bam.log" \
    docker run --rm --name bench_gpu_fq2bam \
        --gpus all \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${BENCH_VOLUME},target=/bench" \
        --mount "type=volume,source=${BENCH_TMP_VOLUME},target=/pbtmp" \
        "$PB_IMAGE" \
        pbrun fq2bam \
            --ref "$REFERENCE" \
            --in-fq "$SUB_R1" "$SUB_R2" "$READ_GROUP" \
            --out-bam /bench/gpu.bam \
            --out-duplicate-metrics /bench/gpu.duplicate_metrics.txt \
            --bwa-options=-Y \
            --low-memory \
            --memory-limit 8 \
            --bwa-normalized-queue-capacity 2 \
            --gpuwrite \
            --monitor-usage \
            --tmp-dir /pbtmp \
            --num-gpus 1
record "A1_gpu_fq2bam" "$TIMED" "align+sort+markdup on GPU"


# -----------------------------------------------------------------------------
# TEST A2 - alignment on CPU: bwa mem | samtools sort, then MarkDuplicates
#
# This is the classic GATK Best Practices path, the one fq2bam replaces. The two
# steps are timed separately because the bottleneck is not the same.
# -----------------------------------------------------------------------------

log "Test A2 - bwa mem + samtools sort (CPU, ${THREADS} threads)"
timed "${LOG_DIR}/A2a_cpu_bwa_sort.log" \
    on_bench bash -c "bwa mem -t ${THREADS} -Y -K 100000000 \
            -R '${READ_GROUP}' \
            '${REFERENCE}' '${SUB_R1}' '${SUB_R2}' \
        | samtools sort -@ 4 -m 1G -T /bench/sorttmp -o /bench/cpu.sorted.bam -"
CPU_ALIGN="$TIMED"
record "A2a_cpu_bwa_sort" "$TIMED" "align+sort on ${THREADS} threads"

log "Test A2 - gatk MarkDuplicates (CPU)"
timed "${LOG_DIR}/A2b_cpu_markdup.log" \
    on_bench gatk --java-options "-Xmx8g" MarkDuplicates \
        -I /bench/cpu.sorted.bam \
        -O /bench/cpu.bam \
        -M /bench/cpu.duplicate_metrics.txt \
        --TMP_DIR /bench
CPU_MARKDUP="$TIMED"
record "A2b_cpu_markdup" "$TIMED" "duplicate marking on CPU"
record "A2_cpu_total" "$((CPU_ALIGN + CPU_MARKDUP))" "CPU path equivalent to fq2bam, total"


# -----------------------------------------------------------------------------
# TEST B1 - variant calling on GPU, full BAM, fixed region
# -----------------------------------------------------------------------------

log "Test B1 - pbrun haplotypecaller (GPU) on ${REGION}"
timed "${LOG_DIR}/B1_gpu_hc.log" \
    docker run --rm --name bench_gpu_hc \
        --gpus all \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${FULL_VOLUME},target=/full,readonly" \
        --mount "type=volume,source=${BENCH_VOLUME},target=/bench" \
        --mount "type=volume,source=${BENCH_TMP_VOLUME},target=/pbtmp" \
        "$PB_IMAGE" \
        pbrun haplotypecaller \
            --ref "$REFERENCE" \
            --in-bam "$FULL_BAM" \
            --out-variants /bench/gpu.region.g.vcf.gz \
            --gvcf \
            --htvc-low-memory \
            -L "$REGION" \
            --tmp-dir /pbtmp \
            --num-gpus 1
record "B1_gpu_haplotypecaller" "$TIMED" "gVCF over ${REGION}, real coverage"


# -----------------------------------------------------------------------------
# TEST B2 - variant calling on CPU: gatk HaplotypeCaller, same BAM, same region,
# same gVCF output.
# -----------------------------------------------------------------------------

log "Test B2 - gatk HaplotypeCaller (CPU) on ${REGION}"
timed "${LOG_DIR}/B2_cpu_hc.log" \
    docker run --rm --name bench_cpu_hc \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${FULL_VOLUME},target=/full,readonly" \
        --mount "type=volume,source=${BENCH_VOLUME},target=/bench" \
        "$TOOLS_IMAGE" \
        gatk --java-options "-Xmx8g" HaplotypeCaller \
            -R "$REFERENCE" \
            -I "$FULL_BAM" \
            -O /bench/cpu.region.g.vcf.gz \
            -ERC GVCF \
            -L "$REGION" \
            --native-pair-hmm-threads "$THREADS"
record "B2_cpu_haplotypecaller" "$TIMED" "gVCF over ${REGION}, real coverage"


# -----------------------------------------------------------------------------
# Concordance: the speed only counts if the two callers say the same thing.
# -----------------------------------------------------------------------------

# No bcftools in the tools container: the gVCFs are read with zcat and awk. In a
# gVCF the block lines have ALT equal to <NON_REF> and are not variants: only
# the ones carrying a real allele are kept.
log "Comparing the two gVCFs"
on_bench bash -c "
    set -e
    for side in gpu cpu; do
        zcat /bench/\${side}.region.g.vcf.gz \
            | awk -F'\t' '!/^#/ && \$5 != \"<NON_REF>\" { print \$1\"\t\"\$2\"\t\"\$4\"\t\"\$5 }' \
            | sort -u > /bench/\${side}.sites.txt
    done
    printf 'gpu_sites\t%s\n' \"\$(wc -l < /bench/gpu.sites.txt)\"
    printf 'cpu_sites\t%s\n' \"\$(wc -l < /bench/cpu.sites.txt)\"
    printf 'shared_sites\t%s\n' \"\$(comm -12 /bench/gpu.sites.txt /bench/cpu.sites.txt | wc -l)\"
    printf 'gpu_only\t%s\n'  \"\$(comm -23 /bench/gpu.sites.txt /bench/cpu.sites.txt | wc -l)\"
    printf 'cpu_only\t%s\n'  \"\$(comm -13 /bench/gpu.sites.txt /bench/cpu.sites.txt | wc -l)\"
" > "${RESULT_DIR}/concordance.tsv" 2> "${LOG_DIR}/C_concordance.log" ||
    log "WARNING: comparison failed, see ${LOG_DIR}/C_concordance.log"


# -----------------------------------------------------------------------------
# Machine-readable summary
# -----------------------------------------------------------------------------

{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date --iso-8601=seconds)"
    printf '  "pairs": %s,\n' "$PAIRS"
    printf '  "region": "%s",\n' "$REGION"
    printf '  "threads": %s,\n' "$THREADS"
    printf '  "gpu": "RTX 3090",\n'
    printf '  "seconds": {\n'
    awk -F'\t' 'NR>1 && $1 != "" {
        rows[n++] = sprintf("    \"%s\": %s", $1, $2)
    }
    END {
        for (i = 0; i < n; i++) printf "%s%s\n", rows[i], (i < n - 1 ? "," : "")
    }' "$RESULT_TSV"
    printf '  }\n'
    printf '}\n'
} > "$RESULT_JSON"

log "Done. Timings in ${RESULT_TSV}"
cat "$RESULT_TSV"
