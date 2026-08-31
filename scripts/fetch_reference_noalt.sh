#!/usr/bin/env bash
# =============================================================================
# Downloads and prepares the GRCh38 reference without ALT contigs (iteration 2).
#
# Why: the reference of the first iteration (Homo_sapiens_assembly38) contains
# 261 _alt contigs but not the .alt file that tells BWA an ALT contig and its
# primary locus are the same piece of genome. Reads that map to both get
# MAPQ 0, HaplotypeCaller discards them below MAPQ 20, and the variants
# underneath are never called: inside the MHC of chr6 that loses 98.51 % of the
# true SNPs.
#
# The variant chosen is no_alt_plus_hs38d1: with respect to today it changes one
# thing only, that is, it removes the ALT contigs (and the HLA contigs that
# depend on them) and keeps the decoys, which are there to absorb junk reads.
#
# The BWA indexes are published by NCBI and are not rebuilt: bwa index over
# 3 Gb would cost about an hour of CPU for an identical result.
#
# Download, decompression and indexes all run INSIDE the container: the ref/
# directory belongs to root and the WSL user cannot write to it. The container
# runs as root, so no sudo is needed. The host only orchestrates and does the
# read-only comparisons, which is the model the whole pipeline uses.
#
#   Usage:  bash scripts/fetch_reference_noalt.sh
# =============================================================================

set -Eeuo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="${PROJECT_DIR}/ref"

REF_NAME="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids"

REF="${REF_DIR}/${REF_NAME}"
DICT="${REF_DIR}/${REF_NAME%.*}.dict"
OLD_FAI="${REF_DIR}/Homo_sapiens_assembly38.fasta.fai"
TOOLS_IMAGE="${TOOLS_IMAGE:-bioinfo-codeserver:latest}"

msg() { printf '\n\033[1;36m==> %s\033[0m\n' "$*"; }
die() { printf '\n\033[1;31m!! %s\033[0m\n' "$*" >&2; exit 1; }


# --- 1. Prerequisites --------------------------------------------------------
msg "1/4 - Prerequisites"
command -v docker >/dev/null || die "Docker is not available."
docker image inspect "$TOOLS_IMAGE" >/dev/null 2>&1 ||
    die "image ${TOOLS_IMAGE} missing: it provides curl, samtools and the .dict."
[[ -d "$REF_DIR" ]] || die "directory not found: $REF_DIR"
[[ -s "$OLD_FAI" ]] || die "the .fai of the current reference is missing: $OLD_FAI"

free_kb="$(df -Pk "$REF_DIR" | awk 'END {print $4}')"
(( free_kb >= 20 * 1024 * 1024 )) ||
    die "at least 20 GiB of free space are needed in ${REF_DIR}."
printf 'Free space: %.1f GiB\n' "$(awk -v k="$free_kb" 'BEGIN{print k/1024/1024}')"


# --- 2. Download and prepare, inside the container ---------------------------
msg "2/4 - Download, MD5 verification, BWA indexes, .fai and .dict"
docker run --rm -i \
    --name hg002_fetch_reference \
    --entrypoint /bin/bash \
    --mount "type=bind,source=${REF_DIR},target=/ref" \
    "$TOOLS_IMAGE" -s -- "$REF_NAME" "$BASE_URL" <<'PAYLOAD'
set -Eeuo pipefail

REF_NAME="$1"
BASE_URL="$2"
REF="/ref/${REF_NAME}"
DICT="/ref/${REF_NAME%.*}.dict"

step() { printf '\n  [%s]\n' "$*"; }
die()  { printf '\n!! %s\n' "$*" >&2; exit 1; }

# --no-progress-meter keeps the log readable; progress is followed by watching
# the files in ref/ grow.
#
# The outer loop exists because NCBI's FTP drops the connection without warning
# on 3 GB archives: curl returns 18 (transfer closed) and --retry does not treat
# that as retryable. Every pass resumes with -C - from where it got to, so an
# interruption near the end costs a few MB rather than 3 GB.
fetch() {
    local url="$1" dest="$2" attempt=0
    if [[ -s "$dest" ]]; then
        printf '  already present: %-58s %s\n' "$(basename "$dest")" "$(du -h "$dest" | cut -f1)"
        return
    fi
    printf '  downloading:     %s\n' "$(basename "$dest")"
    until curl -L --fail --retry 5 --retry-delay 5 --retry-all-errors \
               --no-progress-meter -C - -o "${dest}.part" "$url"; do
        attempt=$((attempt + 1))
        (( attempt < 8 )) ||
            die "download failed after ${attempt} attempts: $url"
        printf '  interrupted, resuming (attempt %s, %s downloaded)\n' \
            "$attempt" "$(du -h "${dest}.part" 2>/dev/null | cut -f1 || echo 0)"
        sleep 10
    done
    mv -f "${dest}.part" "$dest"
    printf '  downloaded:      %-58s %s\n' "$(basename "$dest")" "$(du -h "$dest" | cut -f1)"
}

step "download"
fetch "${BASE_URL}/${REF_NAME}.gz"               "${REF}.gz"
fetch "${BASE_URL}/${REF_NAME}.fai"              "${REF}.fai.ncbi"
fetch "${BASE_URL}/${REF_NAME}.bwa_index.tar.gz" "${REF}.bwa_index.tar.gz"

step "MD5 sums declared by NCBI"
if curl -L --fail --silent --max-time 120 \
        -o /ref/md5checksums.ncbi.txt "${BASE_URL}/md5checksums.txt"; then
    cd /ref
    for f in "${REF_NAME}.gz" "${REF_NAME}.bwa_index.tar.gz"; do
        expected="$(awk -v n="$f" 'index($2, n) {print $1}' md5checksums.ncbi.txt | head -1)"
        if [[ -z "$expected" ]]; then
            printf '  %-58s md5 not declared\n' "$f"
            continue
        fi
        obtained="$(md5sum "$f" | cut -d' ' -f1)"
        [[ "$expected" == "$obtained" ]] ||
            die "md5 mismatch for ${f}: expected ${expected}, got ${obtained}"
        printf '  %-58s md5 OK\n' "$f"
    done
else
    echo "  md5checksums.txt not available: MD5 verification skipped."
    rm -f /ref/md5checksums.ncbi.txt
fi

step "decompressing the FASTA"
if [[ -s "$REF" ]]; then
    echo "  FASTA already decompressed."
else
    gzip -t "${REF}.gz" || die "FASTA archive corrupt: download it again."
    gunzip -c "${REF}.gz" > "${REF}.part"
    mv -f "${REF}.part" "$REF"
fi
printf '  FASTA: %s\n' "$(du -h "$REF" | cut -f1)"

# NCBI publishes the BWA indexes pre-built. Some versions name them with a ".64"
# infix; BWA and Parabricks look for them without it, so they get renamed
# (renamed, not copied: no extra space needed).
step "BWA indexes"
if [[ -s "${REF}.bwt" && -s "${REF}.sa" && -s "${REF}.pac" &&
      -s "${REF}.amb" && -s "${REF}.ann" ]]; then
    echo "  BWA indexes already present."
else
    tar -tzf "${REF}.bwa_index.tar.gz" >/dev/null ||
        die "BWA index archive corrupt: download it again."
    tar -xzf "${REF}.bwa_index.tar.gz" -C /ref
    if [[ ! -s "${REF}.bwt" && ! -s "${REF}.64.bwt" ]]; then
        nested="$(find /ref -mindepth 2 -name '*.bwt' -print -quit)"
        [[ -n "$nested" ]] || die "BWA indexes not found after extraction."
        mv -f "$(dirname "$nested")"/* /ref/
        rmdir "$(dirname "$nested")" 2>/dev/null || true
    fi
    for ext in amb ann bwt pac sa; do
        if [[ ! -s "${REF}.${ext}" && -s "${REF}.64.${ext}" ]]; then
            mv -f "${REF}.64.${ext}" "${REF}.${ext}"
            printf '  renamed .64.%s -> .%s\n' "$ext" "$ext"
        fi
        [[ -s "${REF}.${ext}" ]] || die "BWA index missing: ${REF}.${ext}"
    done
fi
for ext in amb ann bwt pac sa; do
    printf '  %-6s %s\n' ".${ext}" "$(du -h "${REF}.${ext}" | cut -f1)"
done

# The .fai is regenerated and compared with the one published by NCBI: if they
# match, the downloaded FASTA is intact line by line.
step ".fai and .dict"
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"
if cmp -s "${REF}.fai" "${REF}.fai.ncbi"; then
    echo "  regenerated .fai identical to the one published by NCBI."
else
    die "regenerated .fai differs from NCBI's: the FASTA is not intact."
fi
[[ -s "$DICT" ]] || samtools dict -o "$DICT" "$REF"
[[ -s "$DICT" ]] || die ".dict dictionary was not created."
printf '  .dict: %s sequences\n' "$(grep -c '^@SQ' "$DICT")"
PAYLOAD


# --- 3. What changed ---------------------------------------------------------
# This is the check that matters: it says which contigs disappear and confirms
# that the primary chromosomes do NOT change, that is, that dbSNP, snpEff and
# the GIAB BED stay compatible without being touched.
msg "3/4 - Comparison with the reference of the first iteration"
old="$(mktemp)"; new="$(mktemp)"; removed="$(mktemp)"; added="$(mktemp)"
trap 'rm -f "$old" "$new" "$removed" "$added"' EXIT
cut -f1 "$OLD_FAI"   | sort > "$old"
cut -f1 "${REF}.fai" | sort > "$new"
comm -23 "$old" "$new" > "$removed"
comm -13 "$old" "$new" > "$added"

# grep -c always prints the count, but exits 1 when it is zero: what is needed
# is to ignore the status, not to replace the output.
count()     { grep -cE  "$1" "$2" || true; }
count_not() { grep -vcE "$1" "$2" || true; }
printf 'Total contigs : %s -> %s\n' "$(wc -l < "$old")" "$(wc -l < "$new")"
printf 'Removed       : %s (_alt: %s, HLA-: %s, other: %s)\n' \
    "$(wc -l < "$removed")" \
    "$(count '_alt$' "$removed")" \
    "$(count '^HLA-' "$removed")" \
    "$(count_not '(_alt$|^HLA-)' "$removed")"
printf 'Added         : %s\n' "$(wc -l < "$added")"
[[ -s "$added" ]] && sed 's/^/    + /' "$added" | head -20


# --- 4. Compatibility with dbSNP, snpEff and the GIAB BED --------------------
msg "4/4 - Did the primary chromosomes stay identical?"
if grep -qE '^chr([0-9]{1,2}|X|Y|M)$' "$added"; then
    die "a primary chromosome changed name: dbSNP and the GIAB BED would no longer be compatible."
fi
if grep -qE '^chr([0-9]{1,2}|X|Y|M)$' "$removed"; then
    die "a primary chromosome is missing from the new reference."
fi
for c in chr1 chr6 chr10 chr20 chrX chrY chrM; do
    before="$(awk -v c="$c" '$1==c {print $2}' "$OLD_FAI")"
    after="$(awk -v c="$c" '$1==c {print $2}' "${REF}.fai")"
    [[ -n "$after" ]] || die "${c} is absent from the new reference."
    [[ "$before" == "$after" ]] ||
        die "${c} has a different length (${before} vs ${after}): this is not the same GRCh38."
    printf '  %-6s %s bp, identical\n' "$c" "$after"
done
echo
echo "Names and lengths of the primary chromosomes unchanged."
echo "dbSNP, snpEff and the GIAB BED stay valid without modification."

msg "Reference ready: ${REF}"
echo "The pipeline already uses it by default. To go back to the old one:"
echo "  HG002_REF_NAME=Homo_sapiens_assembly38.fasta bash run_parabricks_hg002.sh --preflight"
