#!/usr/bin/env bash
# Run from the repository root in Linux or Ubuntu WSL2.
# Edit the inputs below, then: bash run_parabricks.sh
set -euo pipefail

# Input paths are relative to this directory. Use the read group of your data.
SAMPLE="${SAMPLE:-HG002}"
R1="${R1:-HG002.novaseq.pcr-free.40x.R1.fastq.gz}"
R2="${R2:-HG002.novaseq.pcr-free.40x.R2.fastq.gz}"
REF="ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
KNOWN_SITES="ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
READ_GROUP="${READ_GROUP:-@RG\tID:HV3C3DSXX.2\tPL:ILLUMINA\tPM:NovaSeq6000\tLB:HG002_PCR_FREE\tPU:HV3C3DSXX.2.AGCGATAG+AGGCGAAG\tSM:HG002}"
IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
OUT="${OUT:-output/${SAMPLE}_haplotypecaller}"

# mkdir fails if this run directory already exists. Pick a new OUT to rerun.
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT"
mkdir "$OUT/tmp"

# 1. Align paired reads, sort the BAM and mark duplicates.
docker run --rm --gpus all \
    --volume "$PWD:/input:ro" \
    --volume "$PWD/$OUT:/output" \
    "$IMAGE" pbrun fq2bam \
    --ref "/input/$REF" \
    --in-fq "/input/$R1" "/input/$R2" "$READ_GROUP" \
    --out-bam "/output/$SAMPLE.bam" \
    --out-duplicate-metrics "/output/$SAMPLE.duplicate_metrics.txt" \
    --bwa-options=-Y \
    --low-memory --memory-limit 8 --bwa-normalized-queue-capacity 2 \
    --gpuwrite --monitor-usage --num-gpus 1 --tmp-dir /output/tmp

# 2. Estimate base-quality recalibration in a separate container to release RAM.
docker run --rm --gpus all \
    --volume "$PWD:/input:ro" \
    --volume "$PWD/$OUT:/output" \
    "$IMAGE" pbrun bqsr \
    --ref "/input/$REF" \
    --in-bam "/output/$SAMPLE.bam" \
    --knownSites "/input/$KNOWN_SITES" \
    --out-recal-file "/output/$SAMPLE.recal.txt" \
    --num-gpus 1 --tmp-dir /output/tmp

# 3. Call variants as gVCF, applying recalibration during calling.
# The BAM on disk still has its original base qualities.
docker run --rm --gpus all \
    --volume "$PWD:/input:ro" \
    --volume "$PWD/$OUT:/output" \
    "$IMAGE" pbrun haplotypecaller \
    --ref "/input/$REF" \
    --in-bam "/output/$SAMPLE.bam" \
    --in-recal-file "/output/$SAMPLE.recal.txt" \
    --out-variants "/output/$SAMPLE.g.vcf.gz" \
    --gvcf --htvc-low-memory --num-gpus 1 --tmp-dir /output/tmp

# 4. Convert the gVCF into a single-sample variant VCF.
docker run --rm \
    --volume "$PWD:/input:ro" \
    --volume "$PWD/$OUT:/output" \
    "$IMAGE" pbrun genotypegvcf \
    --ref "/input/$REF" \
    --in-gvcf "/output/$SAMPLE.g.vcf.gz" \
    --out-vcf "/output/$SAMPLE.vcf.gz" \
    --tmp-dir /output/tmp

printf 'Variant calling finished: %s/%s.vcf.gz\n' "$OUT" "$SAMPLE"
