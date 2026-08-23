#!/usr/bin/env bash
# HG002: FASTQ -> aligned BAM -> BQSR -> gVCF -> VCF with NVIDIA Parabricks.
# Large intermediates stay in a Linux Docker volume to avoid WSL/host I/O costs.

set -Eeuo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PB_IMAGE="${PB_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1}"
WORK_VOLUME="${WORK_VOLUME:-hg002_parabricks_work}"
TMP_VOLUME="${TMP_VOLUME:-hg002_parabricks_tmp}"

REF_NAME="${REF_NAME:-GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"
REFERENCE="/project/ref/${REF_NAME}"
KNOWN_SITES="/project/ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
FASTQ_R1="/project/HG002.novaseq.pcr-free.40x.R1.fastq.gz"
FASTQ_R2="/project/HG002.novaseq.pcr-free.40x.R2.fastq.gz"
READ_GROUP='@RG\tID:HV3C3DSXX.2\tPL:ILLUMINA\tPM:NovaSeq6000\tLB:HG002_PCR_FREE\tPU:HV3C3DSXX.2.AGCGATAG+AGGCGAAG\tSM:HG002'

mkdir -p "${PROJECT_DIR}/output"

for file in \
    "${PROJECT_DIR}/ref/${REF_NAME}" \
    "${PROJECT_DIR}/ref/${REF_NAME}.fai" \
    "${PROJECT_DIR}/ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
    "${PROJECT_DIR}/HG002.novaseq.pcr-free.40x.R1.fastq.gz" \
    "${PROJECT_DIR}/HG002.novaseq.pcr-free.40x.R2.fastq.gz"; do
    [[ -s "$file" ]] || { echo "Missing input: $file" >&2; exit 1; }
done

command -v docker >/dev/null || { echo "Docker is required." >&2; exit 1; }
docker info >/dev/null 2>&1 || { echo "Docker is not running." >&2; exit 1; }

# Refuse the exact reference/index mismatch that destroyed MHC recall in run 1.
if grep -q '_alt$' "${PROJECT_DIR}/ref/${REF_NAME}.fai" &&
   [[ ! -s "${PROJECT_DIR}/ref/${REF_NAME}.alt" ]]; then
    echo "Reference contains ALT contigs but has no companion .alt file." >&2
    exit 1
fi

docker volume create "$WORK_VOLUME" >/dev/null
docker volume create "$TMP_VOLUME" >/dev/null

pbrun() {
    docker run --rm --gpus all \
        --mount "type=bind,source=${PROJECT_DIR},target=/project,readonly" \
        --mount "type=volume,source=${WORK_VOLUME},target=/work" \
        --mount "type=volume,source=${TMP_VOLUME},target=/pbtmp" \
        "$PB_IMAGE" pbrun "$@"
}

echo "[1/4] fq2bam"
pbrun fq2bam \
    --ref "$REFERENCE" \
    --in-fq "$FASTQ_R1" "$FASTQ_R2" "$READ_GROUP" \
    --out-bam /work/HG002.bam \
    --out-duplicate-metrics /work/HG002.duplicate_metrics.txt \
    --bwa-options=-Y \
    --low-memory \
    --memory-limit 8 \
    --bwa-normalized-queue-capacity 2 \
    --gpuwrite \
    --tmp-dir /pbtmp \
    --num-gpus 1

echo "[2/4] bqsr"
pbrun bqsr \
    --ref "$REFERENCE" \
    --in-bam /work/HG002.bam \
    --knownSites "$KNOWN_SITES" \
    --out-recal-file /work/HG002.recal.txt \
    --tmp-dir /pbtmp \
    --num-gpus 1

echo "[3/4] haplotypecaller"
pbrun haplotypecaller \
    --ref "$REFERENCE" \
    --in-bam /work/HG002.bam \
    --in-recal-file /work/HG002.recal.txt \
    --out-variants /work/HG002.g.vcf.gz \
    --gvcf \
    --htvc-low-memory \
    --tmp-dir /pbtmp \
    --num-gpus 1

echo "[4/4] genotypegvcf"
pbrun genotypegvcf \
    --ref "$REFERENCE" \
    --in-gvcf /work/HG002.g.vcf.gz \
    --out-vcf /work/HG002.vcf.gz \
    --tmp-dir /pbtmp

# Only the compact final call set is copied out; BAM and gVCF stay in the volume.
docker run --rm \
    --mount "type=volume,source=${WORK_VOLUME},target=/work,readonly" \
    --mount "type=bind,source=${PROJECT_DIR}/output,target=/output" \
    "$PB_IMAGE" bash -c \
    'cp /work/HG002.vcf.gz /output/ && cp /work/HG002.vcf.gz.tbi /output/'

docker volume rm "$TMP_VOLUME" >/dev/null
echo "Done: ${PROJECT_DIR}/output/HG002.vcf.gz"
