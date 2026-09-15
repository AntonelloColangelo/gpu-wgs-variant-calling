#!/usr/bin/env bash
# Run from the repository root after producing the BAM with run_parabricks.sh.
set -euo pipefail

SAMPLE="${SAMPLE:-HG003}"
BAM="${BAM:-output/${SAMPLE}_haplotypecaller/${SAMPLE}.bam}"
REF="ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
OUT="${OUT:-output/${SAMPLE}_deepvariant}"
IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"

# Inputs and OUT are relative to the repository root.
# Use the same BAM as HaplotypeCaller, with original base qualities.
# Use a fresh directory to avoid overwriting another run.
mkdir -p "$(dirname "$OUT")"
mkdir "$OUT"
mkdir "$OUT/tmp"

docker run --rm --gpus all \
    --volume "$PWD:/input:ro" \
    --volume "$PWD/$OUT:/output" \
    "$IMAGE" pbrun deepvariant \
    --ref "/input/$REF" \
    --in-bam "/input/$BAM" \
    --out-variants "/output/$SAMPLE.deepvariant.vcf.gz" \
    --num-gpus 1 \
    --tmp-dir /output/tmp

printf 'Variant calling finished: %s/%s.deepvariant.vcf.gz\n' "$OUT" "$SAMPLE"
