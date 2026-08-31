#!/usr/bin/env bash
# =============================================================================
# HG002 germline pipeline
#
# It reads top to bottom:
#   1. check the FASTQ files
#   2. build and check the BAM
#   3. compute BQSR and call the variants
#   4. filter and annotate the VCF
#   5. produce QC, export and report
# =============================================================================

set -Eeuo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${PROJECT_DIR}/scripts/pipeline_functions.sh"


# -----------------------------------------------------------------------------
# 0. Run mode
# -----------------------------------------------------------------------------

MODE="full"
PREFLIGHT_ONLY="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --smoke)     MODE="smoke" ;;
        --preflight) PREFLIGHT_ONLY="true" ;;
        --help|-h)
            echo "Usage:"
            echo "  ./run_parabricks_hg002.sh              # full 40x analysis"
            echo "  ./run_parabricks_hg002.sh --smoke      # small end-to-end trial"
            echo "  ./run_parabricks_hg002.sh --preflight  # checks only"
            exit 0
            ;;
        *) die "unknown option: $1" ;;
    esac
    shift
done


# -----------------------------------------------------------------------------
# 1. Project data
# -----------------------------------------------------------------------------

SAMPLE="HG002"

# Reference of iteration 2: GRCh38 without ALT contigs, decoys kept.
# The one used in iteration 1 carried 261 _alt contigs without the .alt file
# beside them, and that cost 98.51 % of the true SNPs in the MHC of chr6.
# Download it with: bash scripts/fetch_reference_noalt.sh
#
# To go back to the old reference without editing this script:
#   HG002_REF_NAME=Homo_sapiens_assembly38.fasta \
#   HG002_ALLOW_ALT_WITHOUT_ALT_FILE=1 bash run_parabricks_hg002.sh --preflight
REF_NAME="${HG002_REF_NAME:-GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"
KNOWN_SITES_NAME="Homo_sapiens_assembly38.dbsnp138.vcf.gz"

# Two names for the same file: the host path is what the checks and the
# provenance key use, the container path is what pbrun and gatk get.
REFERENCE_HOST="${PROJECT_DIR}/ref/${REF_NAME}"
REFERENCE_DICT_HOST="${PROJECT_DIR}/ref/${REF_NAME%.*}.dict"
KNOWN_SITES_HOST="${PROJECT_DIR}/ref/${KNOWN_SITES_NAME}"
REFERENCE="/project/ref/${REF_NAME}"
REFERENCE_DICT="/project/ref/${REF_NAME%.*}.dict"
KNOWN_SITES="/project/ref/${KNOWN_SITES_NAME}"

# Read group taken from the NovaSeq header:
# flowcell HV3C3DSXX, lane 2, indexes AGCGATAG+AGGCGAAG.
READ_GROUP_PU="HV3C3DSXX.2.AGCGATAG+AGGCGAAG"
READ_GROUP='@RG\tID:HV3C3DSXX.2\tPL:ILLUMINA\tPM:NovaSeq6000\tLB:HG002_PCR_FREE\tPU:HV3C3DSXX.2.AGCGATAG+AGGCGAAG\tSM:HG002'

PB_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
TOOLS_IMAGE="bioinfo-codeserver:latest"

if [[ "$MODE" == "smoke" ]]; then
    BASE_PREFIX="HG002_NovaSeq_smoke_test"
    FASTQ_R1_HOST="${PROJECT_DIR}/smoke/HG002_NovaSeq_smoke_R1.fastq.gz"
    FASTQ_R2_HOST="${PROJECT_DIR}/smoke/HG002_NovaSeq_smoke_R2.fastq.gz"
    FASTQ_R1_CONTAINER="/project/smoke/HG002_NovaSeq_smoke_R1.fastq.gz"
    FASTQ_R2_CONTAINER="/project/smoke/HG002_NovaSeq_smoke_R2.fastq.gz"
    EXPECTED_R1_BYTES="75559067"
    EXPECTED_R2_BYTES="77823079"
    EXPECTED_PAIRS="1000000"
    BASE_WORK_VOLUME="hg002_smoke_work"
    BASE_TMP_VOLUME="hg002_smoke_tmp"
    INTERVAL_ARGS=(-L "chr20:1-1000000")
else
    BASE_PREFIX="HG002_NovaSeq_40x"
    FASTQ_R1_HOST="${PROJECT_DIR}/HG002.novaseq.pcr-free.40x.R1.fastq.gz"
    FASTQ_R2_HOST="${PROJECT_DIR}/HG002.novaseq.pcr-free.40x.R2.fastq.gz"
    FASTQ_R1_CONTAINER="/project/HG002.novaseq.pcr-free.40x.R1.fastq.gz"
    FASTQ_R2_CONTAINER="/project/HG002.novaseq.pcr-free.40x.R2.fastq.gz"
    EXPECTED_R1_BYTES="34156056301"
    EXPECTED_R2_BYTES="35449847992"
    EXPECTED_PAIRS="474384500"
    BASE_WORK_VOLUME="hg002_work"
    BASE_TMP_VOLUME="hg002_tmp"
    INTERVAL_ARGS=()
fi


# -----------------------------------------------------------------------------
# 1b. Provenance key
#
# Prefix and volumes are not hand-written constants: they derive from everything
# that determines the result (container image, reference, FASTQ files,
# known-sites, read group, intervals). Changing one of those ingredients
# produces new names, so the checkpoints of the previous run are neither reused
# by mistake nor overwritten.
#
# Without this, changing the reference and re-launching would find the old BAM
# again, and it would pass validation: the pipeline would report PASS while
# returning exactly the results it was supposed to replace.
# -----------------------------------------------------------------------------

require_runtime
RUN_KEY="$(compute_run_key)"
PREFIX="${BASE_PREFIX}_${RUN_KEY}"
WORK_VOLUME="${BASE_WORK_VOLUME}_${RUN_KEY}"
TMP_VOLUME="${BASE_TMP_VOLUME}_${RUN_KEY}"


# -----------------------------------------------------------------------------
# 2. Logging, and protection against two simultaneous starts
# -----------------------------------------------------------------------------

OUTPUT_DIR="${PROJECT_DIR}/output"
REPORT_DIR="${PROJECT_DIR}/reports"
LOG_DIR="${PROJECT_DIR}/logs"
RUN_ID="$(date '+%Y%m%d_%H%M%S')_${MODE}"
RUN_DIR="${LOG_DIR}/runs/${RUN_ID}"
TIMING_FILE="${RUN_DIR}/${PREFIX}.step_times.tsv"
STATE_FILE="${LOG_DIR}/current_step.env"

mkdir -p "$OUTPUT_DIR" "$REPORT_DIR" "$RUN_DIR"

exec 9>"${LOG_DIR}/pipeline.lock"
flock -n 9 || die "another HG002 pipeline is already running."

cat > "${LOG_DIR}/current_run.env.tmp" <<EOF
RUN_ID=${RUN_ID}
MODE=${MODE}
PREFIX=${PREFIX}
WORK_VOLUME=${WORK_VOLUME}
RUN_KEY=${RUN_KEY}
STARTED=$(date --iso-8601=seconds)
EOF
mv -f "${LOG_DIR}/current_run.env.tmp" "${LOG_DIR}/current_run.env"

# The manifest lists the ingredients of the key in plain text. It is what lets
# you verify, months later, where a result came from: a hash on its own tells
# you that something changed, not what.
write_run_manifest "${RUN_DIR}/run_manifest.json"

printf 'Step\tStarted\tEnded\tSeconds\tStatus\n' > "$TIMING_FILE"
exec 3>&1 4>&2
exec > >(tee -a "${RUN_DIR}/pipeline.log" >&3) 2>&1
TEE_PID="$!"

on_pipeline_exit() {
    local rc="$?"
    if [[ "$rc" -ne 0 ]]; then
        pipeline_failed "$rc" "${BASH_LINENO[0]:-?}"
    fi
    # Closes the stream to tee and waits for the log to be written in full.
    exec 1>&3 2>&4
    wait "$TEE_PID" 2>/dev/null || true
    return "$rc"
}
trap on_pipeline_exit EXIT

# The FASTQ integrity check depends only on the FASTQ files, not on the
# reference: a checkpoint from a previous run stays valid and redoing it would
# cost hours for nothing. That is why the list also includes the prefixes from
# before the provenance key existed.
VALIDATION_CANDIDATES=(
    "${OUTPUT_DIR}/${PREFIX}.fastq_validation.json"
    "${OUTPUT_DIR}/${BASE_PREFIX}.fastq_validation.json"
    "${OUTPUT_DIR}/HG002_NovaSeq_${MODE}.fastq_validation.json"
    "${LOG_DIR}/previous_attempts/HG002_NovaSeq_${MODE}.fastq_validation.json"
    "${LOG_DIR}/previous_attempts/HG002_NovaSeq_40x.fastq_validation.json"
    "${LOG_DIR}/previous_attempts/HG002_NovaSeq_smoke.fastq_validation.json"
)

echo "============================================================"
echo " HG002 pipeline | mode: ${MODE}"
echo " Run: ${RUN_ID}"
echo " Reference: ${REF_NAME}"
echo " Provenance key: ${RUN_KEY}"
echo " Prefix: ${PREFIX}"
echo "============================================================"


# =============================================================================
# STAGE 1 - Preflight checks and FASTQ integrity
# =============================================================================

preflight
ensure_fastq_integrity

if [[ "$PREFLIGHT_ONLY" == "true" ]]; then
    write_state "preflight" "Preflight checks" "PASS" "$(date --iso-8601=seconds)"
    echo
    echo "PREFLIGHT PASSED. No BAM or VCF was computed."
    exit 0
fi


# =============================================================================
# STAGE 2 - FASTQ -> BAM, aligned, sorted and with duplicates marked
# =============================================================================

if valid_bam; then
    skip_step "2_fq2bam" "2/9 - fq2bam"
    skip_step "2b_bam_check" "BAM integrity check"
else
    if valid_partial_bam &&
       work_file_exists "/work/${PREFIX}.duplicate_metrics.partial.txt"; then
        skip_step "2_fq2bam" "2/9 - fq2bam (partial output already complete)"
    else
        work_remove \
            "/work/${PREFIX}.partial.bam" \
            "/work/${PREFIX}.partial.bam.bai" \
            "/work/${PREFIX}.duplicate_metrics.partial.txt"

        run_step "2_fq2bam" "2/9 - fq2bam: alignment, sort and duplicates" \
            run_parabricks "hg002_2_fq2bam" \
            pbrun fq2bam \
            --ref "$REFERENCE" \
            --in-fq "$FASTQ_R1_CONTAINER" "$FASTQ_R2_CONTAINER" "$READ_GROUP" \
            --out-bam "/work/${PREFIX}.partial.bam" \
            --out-duplicate-metrics "/work/${PREFIX}.duplicate_metrics.partial.txt" \
            --bwa-options=-Y \
            --low-memory \
            --memory-limit 8 \
            --bwa-normalized-queue-capacity 2 \
            --gpuwrite \
            --monitor-usage \
            --tmp-dir /pbtmp \
            --num-gpus 1
    fi

    run_step "2b_bam_check" "BAM check: structure, index and read group" \
        finalize_bam
    valid_bam || die "the final BAM does not pass the check."
fi


# =============================================================================
# STAGE 3 - BQSR on its own, then variant calling
# =============================================================================

# BQSR is separate from fq2bam: that way the RAM is released between the two
# processes.
if valid_recal; then
    skip_step "3_bqsr" "3/9 - Base Quality Score Recalibration"
    skip_step "3b_bqsr_check" "BQSR report check"
else
    if tools_quiet bash /project/scripts/postprocess.sh check-recal \
        "/work/${PREFIX}.recal.partial.txt"; then
        skip_step "3_bqsr" "3/9 - BQSR (partial report already complete)"
    else
        work_remove "/work/${PREFIX}.recal.partial.txt"
        run_step "3_bqsr" "3/9 - BQSR: quality recalibration model" \
            run_parabricks "hg002_3_bqsr" \
            pbrun bqsr \
            --ref "$REFERENCE" \
            --in-bam "/work/${PREFIX}.bam" \
            --knownSites "$KNOWN_SITES" \
            --out-recal-file "/work/${PREFIX}.recal.partial.txt" \
            "${INTERVAL_ARGS[@]}" \
            --tmp-dir /pbtmp \
            --num-gpus 1
    fi
    run_step "3b_bqsr_check" "BQSR report check" finalize_recal
    valid_recal || die "the final BQSR report is not valid."
fi

# HaplotypeCaller produces a gVCF and applies BQSR on the fly with
# --in-recal-file.
GVCF="/work/${PREFIX}.g.vcf.gz"
GVCF_PARTIAL="/work/${PREFIX}.partial.g.vcf.gz"
if valid_vcf "$GVCF"; then
    skip_step "4_haplotypecaller" "4/9 - HaplotypeCaller (gVCF)"
    skip_step "4b_gvcf_check" "gVCF integrity check"
else
    if valid_vcf "$GVCF_PARTIAL"; then
        skip_step "4_haplotypecaller" "4/9 - HaplotypeCaller (partial gVCF valid)"
    else
        work_remove "$GVCF_PARTIAL" "${GVCF_PARTIAL}.tbi"
        run_step "4_haplotypecaller" "4/9 - HaplotypeCaller: producing the gVCF" \
            run_parabricks "hg002_4_haplotypecaller" \
            pbrun haplotypecaller \
            --ref "$REFERENCE" \
            --in-bam "/work/${PREFIX}.bam" \
            --in-recal-file "/work/${PREFIX}.recal.txt" \
            --out-variants "$GVCF_PARTIAL" \
            --gvcf \
            --htvc-low-memory \
            "${INTERVAL_ARGS[@]}" \
            --tmp-dir /pbtmp \
            --num-gpus 1
    fi
    run_step "4b_gvcf_check" "gVCF check and finalisation" \
        finalize_vcf "$GVCF_PARTIAL" "$GVCF" "hg002_4b_gvcf_check"
    valid_vcf "$GVCF" || die "the final gVCF is not valid."
fi

# GenotypeGVCF turns the gVCF into an ordinary variant VCF.
VCF="/work/${PREFIX}.vcf.gz"
VCF_PARTIAL="/work/${PREFIX}.partial.vcf.gz"
if valid_vcf "$VCF"; then
    skip_step "5_genotypegvcf" "5/9 - GenotypeGVCF"
    skip_step "5b_vcf_check" "VCF integrity check"
else
    if valid_vcf "$VCF_PARTIAL"; then
        skip_step "5_genotypegvcf" "5/9 - GenotypeGVCF (partial VCF valid)"
    else
        work_remove "$VCF_PARTIAL" "${VCF_PARTIAL}.tbi"
        run_step "5_genotypegvcf" "5/9 - GenotypeGVCF: gVCF -> VCF" \
            run_parabricks "hg002_5_genotypegvcf" \
            pbrun genotypegvcf \
            --ref "$REFERENCE" \
            --in-gvcf "$GVCF" \
            --out-vcf "$VCF_PARTIAL" \
            --tmp-dir /pbtmp
    fi
    run_step "5b_vcf_check" "VCF check and finalisation" \
        finalize_vcf "$VCF_PARTIAL" "$VCF" "hg002_5b_vcf_check"
    valid_vcf "$VCF" || die "the final VCF is not valid."
fi


# =============================================================================
# STAGE 4 - Hard filtering and snpEff annotation
# =============================================================================

if valid_hardfilter; then
    skip_step "6_hardfilter" "6/9 - Hard filtering"
else
    run_step "6_hardfilter" "6/9 - Hard filtering without removing records" \
        run_tools "hg002_6_hardfilter" \
        bash /project/scripts/postprocess.sh hardfilter "$PREFIX" "$SAMPLE"
    valid_hardfilter || die "the hard-filtered VCF is not valid."
fi

if valid_annotation; then
    skip_step "7_annotation" "7/9 - snpEff annotation"
else
    run_step "7_annotation" "7/9 - snpEff annotation and PASS+HIGH table" \
        run_tools "hg002_7_annotation" \
        bash /project/scripts/postprocess.sh annotate "$PREFIX" "$SAMPLE"
    valid_annotation || die "the annotation outputs are not valid."
fi


# =============================================================================
# STAGE 5 - Full QC, export to Windows and HTML report
# =============================================================================

if valid_qc; then
    skip_step "8_qc" "8/9 - Full QC"
else
    run_step "8_qc" "8/9 - QC: alignment, coverage and WGS metrics" \
        run_tools "hg002_8_qc" \
        bash /project/scripts/postprocess.sh qc "$PREFIX" "$MODE"
    valid_qc || die "the required QC metrics are incomplete."
fi

# The large files are copied to C: only now, after every check. If the copy
# failed, the valid originals would stay in the Linux volume.
run_step "8b_export" "Exporting the already validated results to Windows" \
    export_results

run_step "9_report" "9/9 - HTML report and JSON summary" \
    generate_report

[[ -s "${REPORT_DIR}/${PREFIX}_Report.html" ]] ||
    die "the HTML report was not created."
[[ -s "${REPORT_DIR}/${PREFIX}_summary.json" ]] ||
    die "the JSON summary was not created."

# The temporary volume belongs to this pipeline and holds no results.
docker volume rm "$TMP_VOLUME" >/dev/null 2>&1 || true

write_state "complete" "Pipeline complete" "PASS" "$(date --iso-8601=seconds)"

echo
echo "============================================================"
echo " PIPELINE COMPLETE AND VALIDATED"
echo "============================================================"
echo " Results: ${OUTPUT_DIR}"
echo " Report:  ${REPORT_DIR}/${PREFIX}_Report.html"
echo " Logs:    ${RUN_DIR}"
echo " Recovery volume: ${WORK_VOLUME}"
echo
echo " The GIAB benchmark was not run: it stays separate."
echo "============================================================"
