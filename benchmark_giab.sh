#!/usr/bin/env bash
# =============================================================================
# GIAB benchmark: how accurate are the variants this pipeline calls?
#
# Compares the VCF the pipeline produced with the NIST/GIAB truth set for HG002 or HG003
# (v4.2.1, GRCh38) using hap.py with the vcfeval engine.
#
# The comparison is restricted to the high-confidence BED: outside those regions
# not even GIAB knows what the right answer is, so counting them as errors would
# skew the result.
#
# It stays deliberately separate from the main pipeline: it runs once, after a
# run has finished, and touches no result.
#
#   Usage:  bash benchmark_giab.sh
# =============================================================================

set -Eeuo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The reference has to be the one the VCF was called against: vcfeval normalises
# variants by re-reading the bases, and with the wrong reference the counts are
# meaningless.
REF_NAME="${HG002_REF_NAME:-GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"
REFERENCE="${PROJECT_DIR}/ref/${REF_NAME}"
REFERENCE_CONTAINER="/data/ref/${REF_NAME}"

# Select the sample, query and result label explicitly.
SAMPLE="${SAMPLE:-HG002}"
case "$SAMPLE" in
    HG002) GIAB_SAMPLE="HG002_NA24385_son" ;;
    HG003) GIAB_SAMPLE="HG003_NA24149_father" ;;
    *) printf 'SAMPLE must be HG002 or HG003.\n' >&2; exit 1 ;;
esac
QUERY_VCF="${QUERY_VCF:-output/${SAMPLE}_haplotypecaller/${SAMPLE}.vcf.gz}"
[[ "$QUERY_VCF" = /* ]] || QUERY_VCF="${PROJECT_DIR}/${QUERY_VCF}"
QUERY_REL="$(realpath --relative-to="$PROJECT_DIR" "$QUERY_VCF")"
[[ "$QUERY_REL" != ../* ]] || { echo 'QUERY_VCF must be inside the repository.' >&2; exit 1; }
BENCH_LABEL="${BENCH_LABEL:-${SAMPLE}_haplotypecaller}"
[[ "$BENCH_LABEL" =~ ^[A-Za-z0-9_.-]+$ && "$BENCH_LABEL" != . && "$BENCH_LABEL" != .. ]] || {
    echo 'BENCH_LABEL must be a simple directory name.' >&2; exit 1;
}
GIAB_DIR="${PROJECT_DIR}/ref/giab"
SDF_DIR="${GIAB_DIR}/${REF_NAME%.*}.sdf"
OUT_DIR="${PROJECT_DIR}/reports/giab/${BENCH_LABEL}"
LOG_FILE="${OUT_DIR}/happy.log"
HAPPY_IMG="jmcdani20/hap.py:v0.3.12"
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/${GIAB_SAMPLE}/NISTv4.2.1/GRCh38"
TRUTH_VCF="${GIAB_DIR}/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="${GIAB_DIR}/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

DOCKER="${DOCKER:-docker}"

msg() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31m!! %s\033[0m\n' "$*" >&2; exit 1; }

mkdir -p "$GIAB_DIR" "$(dirname "$OUT_DIR")"
mkdir "$OUT_DIR" || die "Choose a new BENCH_LABEL: ${OUT_DIR} already exists."

# --- 1. Prerequisites --------------------------------------------------------
msg "1/5 - Checking prerequisites"
command -v "$DOCKER" >/dev/null || die "Docker is not available."
[[ -s "$QUERY_VCF" ]] || die "VCF to evaluate not found: $QUERY_VCF
   Run bash run_parabricks.sh first, or set QUERY_VCF to an existing VCF."
[[ -s "$REFERENCE" ]] || die "Reference not found: $REFERENCE"

if ! $DOCKER image inspect "$HAPPY_IMG" >/dev/null 2>&1; then
    echo "Image $HAPPY_IMG missing: pulling it..."
    $DOCKER pull "$HAPPY_IMG"
fi
echo "VCF to evaluate: ${QUERY_REL} ($(du -h "$QUERY_VCF" | cut -f1))"
echo "Reference:       ${REF_NAME}"

# --- 2. Truth set ------------------------------------------------------------
msg "2/5 - GIAB v4.2.1 truth set (about 250 MB, first time only)"
fetch() {
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then
        echo "Already present: $(basename "$dest")"
        return
    fi
    echo "Downloading $(basename "$dest") ..."
    curl -L --fail --retry 3 -C - -o "$dest" "$url" || die "Download failed: $url"
}
fetch "${GIAB_BASE}/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"                  "$TRUTH_VCF"
fetch "${GIAB_BASE}/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"              "${TRUTH_VCF}.tbi"
fetch "${GIAB_BASE}/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"      "$TRUTH_BED"

# --- 3. RTG index of the reference (SDF) -------------------------------------
# vcfeval works on its own index of the reference. Building it once makes later
# runs much faster. If the build fails it is not serious: hap.py creates a
# temporary one for itself on every pass.
msg "3/5 - RTG index of the reference (one-off, a few minutes)"
if [[ -d "$SDF_DIR" && -n "$(ls -A "$SDF_DIR" 2>/dev/null)" ]]; then
    echo "SDF already present, skipping."
    USE_SDF=1
else
    rm -rf "$SDF_DIR"
    if $DOCKER run --rm \
            --volume "${PROJECT_DIR}":/data \
            --env "REF_CONTAINER=${REFERENCE_CONTAINER}" \
            --env "SDF_CONTAINER=/data/ref/giab/$(basename "$SDF_DIR")" \
            "$HAPPY_IMG" \
            bash -lc 'RTG="$(command -v rtg || true)"
                      [[ -n "$RTG" ]] || RTG=/opt/hap.py/libexec/rtg-tools-install/rtg
                      [[ -x "$RTG" ]] || exit 42
                      "$RTG" format -o "$SDF_CONTAINER" "$REF_CONTAINER"'; then
        echo "SDF created."
        USE_SDF=1
    else
        echo "Could not create the SDF: hap.py will build one for itself."
        rm -rf "$SDF_DIR"
        USE_SDF=0
    fi
fi

# --- 4. hap.py ---------------------------------------------------------------
# A single run is enough: hap.py reports results both for the whole call set
# (the ALL row) and for the PASS variants alone (the PASS row).
msg "4/5 - Comparison with hap.py (vcfeval engine)"

CORES="$(nproc)"
THREADS=$(( CORES < 8 ? CORES : 8 ))   # with 32 GB of RAM, better not to overdo it
echo "Using ${THREADS} threads out of ${CORES} available."

SDF_OPTIONS=()
[[ "$USE_SDF" -eq 1 ]] &&
    SDF_OPTIONS=(--engine-vcfeval-template "/data/ref/giab/$(basename "$SDF_DIR")")

$DOCKER run --rm \
    --volume "${PROJECT_DIR}":/data \
    --env "HGREF=${REFERENCE_CONTAINER}" \
    --env RTG_MEM=12G \
    "$HAPPY_IMG" \
    /opt/hap.py/bin/hap.py \
        "/data/ref/giab/$(basename "$TRUTH_VCF")" \
        "/data/${QUERY_REL}" \
        -f "/data/ref/giab/$(basename "$TRUTH_BED")" \
        -r "$REFERENCE_CONTAINER" \
        -o "/data/reports/giab/${BENCH_LABEL}/happy" \
        --engine=vcfeval \
        "${SDF_OPTIONS[@]}" \
        --threads "$THREADS" \
    2>&1 | tee "$LOG_FILE"

[[ -s "${OUT_DIR}/happy.summary.csv" ]] || die "hap.py produced no summary. See $LOG_FILE"

# --- 5. Readable table -------------------------------------------------------
msg "5/5 - Results"
python3 - "$OUT_DIR" "$SAMPLE" <<'PY'
import csv, json, sys, os

out_dir = sys.argv[1]
rows = list(csv.DictReader(open(os.path.join(out_dir, "happy.summary.csv"))))

def num(value, digits=4):
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return None

header = f"{'Type':<7}{'Filter':<8}{'Recall':>9}{'Precision':>11}{'F1':>9}{'TP':>10}{'FN':>9}{'FP':>9}"
print(header)
print("-" * len(header))

results = {}
for row in rows:
    kind, kept = row["Type"], row["Filter"]
    rec, prec, f1 = num(row["METRIC.Recall"]), num(row["METRIC.Precision"]), num(row["METRIC.F1_Score"])
    print(f"{kind:<7}{kept:<8}"
          f"{rec if rec is not None else 'n/a':>9}"
          f"{prec if prec is not None else 'n/a':>11}"
          f"{f1 if f1 is not None else 'n/a':>9}"
          f"{row['TRUTH.TP']:>10}{row['TRUTH.FN']:>9}{row['QUERY.FP']:>9}")
    results[f"{kind}_{kept}"] = {
        "recall": rec, "precision": prec, "f1": f1,
        "truth_total": num(row["TRUTH.TOTAL"], 0), "tp": num(row["TRUTH.TP"], 0),
        "fn": num(row["TRUTH.FN"], 0), "fp": num(row["QUERY.FP"], 0),
    }

with open(os.path.join(out_dir, "giab_benchmark.json"), "w", encoding="utf-8") as f:
    json.dump({"truth_set": f"{sys.argv[2]} GIAB v4.2.1 GRCh38 (chr1-22, high-confidence BED)",
               "engine": "vcfeval", "results": results}, f, indent=2)
print(f"\nJSON summary: {os.path.join(out_dir, 'giab_benchmark.json')}")
PY

msg "Benchmark complete. Details in ${OUT_DIR}"
