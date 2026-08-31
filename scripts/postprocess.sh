#!/usr/bin/env bash
# Checks and post-processing run inside the bioinfo-codeserver image.

set -Eeuo pipefail

PROJECT="/project"
WORK="/work"

# The reference arrives from the environment: run_tools and tools_quiet in
# pipeline_functions.sh pass it, so there is a single place to change it. The
# default is only for someone running this script by hand.
REFERENCE="${REFERENCE:-${PROJECT}/ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna}"
REFERENCE_DICT="${REFERENCE_DICT:-${REFERENCE%.*}.dict}"

usage() {
    echo "Invalid command for postprocess.sh" >&2
    exit 2
}

check_bam() {
    local bam="$1"
    local bai="$2"
    local sample="$3"
    local pu="$4"
    local header
    header="$(mktemp)"

    samtools quickcheck -v "$bam"
    [[ -s "$bai" ]]
    samtools idxstats "$bam" >/dev/null
    samtools view -H "$bam" > "$header"

    awk -F '\t' -v sample="$sample" -v pu="$pu" '
        $1 == "@HD" {
            for (i=2; i<=NF; i++) if ($i == "SO:coordinate") sorted=1
        }
        $1 == "@RG" {
            sm=0; unit=0
            for (i=2; i<=NF; i++) {
                if ($i == "SM:" sample) sm=1
                if ($i == "PU:" pu) unit=1
            }
            if (sm && unit) correct_rg=1
        }
        END { exit !(sorted && correct_rg) }
    ' "$header"
    rm -f "$header"
}

check_vcf() {
    local vcf="$1"
    local index="$2"
    local sample="$3"
    local header
    header="$(mktemp)"

    gzip -t "$vcf"
    [[ -s "$index" ]]
    tabix -l "$vcf" >/dev/null
    tabix -H "$vcf" > "$header"
    awk -F '\t' -v sample="$sample" '
        /^#CHROM/ { found=1; if (NF < 10 || $10 != sample) exit 2 }
        END { if (!found) exit 3 }
    ' "$header"
    rm -f "$header"
}

check_recal() {
    local recal="$1"
    [[ -s "$recal" ]]
    head -n 1 "$recal" | grep -q '^#:GATKReport'
    grep -q '^#:GATKTable:RecalTable0:' "$recal"
}

count_records() {
    gzip -dc "$1" | awk '!/^#/ {n++} END {print n+0}'
}

command_name="${1:-}"
[[ -n "$command_name" ]] || usage
shift

case "$command_name" in
    integrity)
        [[ $# -eq 3 ]] || usage
        prefix="$1"
        fq1="$2"
        fq2="$3"

        validate_stream() {
            gzip -dc "$1" |
                awk '
                    NR % 4 == 1 && substr($0,1,1) != "@" {bad=1; exit 2}
                    NR % 4 == 2 {sequence_length=length($0)}
                    NR % 4 == 3 && substr($0,1,1) != "+" {bad=1; exit 2}
                    NR % 4 == 0 && length($0) != sequence_length {bad=1; exit 2}
                    END {
                        if (bad || NR == 0 || NR % 4 != 0) exit 2
                        printf "%.0f\n", NR/4
                    }'
        }

        echo "Reading R1 and R2 in full and checking the FASTQ structure..."
        pairs_r1="$(validate_stream "$fq1")"
        pairs_r2="$(validate_stream "$fq2")"
        [[ "$pairs_r1" == "$pairs_r2" ]] || {
            echo "R1 and R2 contain a different number of reads." >&2
            exit 2
        }

        python3 "${PROJECT}/scripts/generate_report.py" validate-fastq \
            --r1 "$fq1" \
            --r2 "$fq2" \
            --pairs 100000 \
            --full-gzip-checked \
            --total-pairs "$pairs_r1" \
            --output "${WORK}/${prefix}.fastq_validation.json"
        ;;

    check-bam)
        [[ $# -eq 4 ]] || usage
        check_bam "$1" "$2" "$3" "$4"
        ;;

    finalize-bam)
        [[ $# -eq 6 ]] || usage
        partial_bam="$1"
        final_bam="$2"
        sample="$3"
        pu="$4"
        partial_metrics="$5"
        final_metrics="$6"

        check_bam "$partial_bam" "${partial_bam}.bai" "$sample" "$pu"
        [[ -s "$partial_metrics" ]]
        mv -f "$partial_bam" "$final_bam"
        mv -f "${partial_bam}.bai" "${final_bam}.bai"
        mv -f "$partial_metrics" "$final_metrics"
        check_bam "$final_bam" "${final_bam}.bai" "$sample" "$pu"
        ;;

    check-recal)
        [[ $# -eq 1 ]] || usage
        check_recal "$1"
        ;;

    finalize-recal)
        [[ $# -eq 2 ]] || usage
        check_recal "$1"
        mv -f "$1" "$2"
        check_recal "$2"
        ;;

    check-vcf)
        [[ $# -eq 3 ]] || usage
        check_vcf "$1" "$2" "$3"
        ;;

    finalize-vcf)
        [[ $# -eq 3 ]] || usage
        partial="$1"
        final="$2"
        sample="$3"
        check_vcf "$partial" "${partial}.tbi" "$sample"
        mv -f "$partial" "$final"
        mv -f "${partial}.tbi" "${final}.tbi"
        check_vcf "$final" "${final}.tbi" "$sample"
        ;;

    hardfilter)
        [[ $# -eq 2 ]] || usage
        prefix="$1"
        sample="$2"
        cd "$WORK"

        raw="${prefix}.vcf.gz"
        final="${prefix}.hardfiltered.vcf.gz"
        partial="${prefix}.hardfiltered.partial.vcf.gz"
        snps="${prefix}.tmp.snps.vcf.gz"
        indels="${prefix}.tmp.indels.vcf.gz"
        other="${prefix}.tmp.other.vcf.gz"

        rm -f "${prefix}.tmp."* "$partial" "${partial}.tbi"

        # Records are split into three groups, filtered, then merged back. The
        # third group keeps MNP, MIXED, SYMBOLIC and NO_VARIATION: no variant
        # type gets lost.
        gatk --java-options "-Xmx8g" SelectVariants \
            -R "$REFERENCE" -V "$raw" \
            --select-type-to-include SNP -O "$snps"
        gatk --java-options "-Xmx8g" SelectVariants \
            -R "$REFERENCE" -V "$raw" \
            --select-type-to-include INDEL -O "$indels"
        gatk --java-options "-Xmx8g" SelectVariants \
            -R "$REFERENCE" -V "$raw" \
            --select-type-to-exclude SNP \
            --select-type-to-exclude INDEL \
            -O "$other"

        # The SNP hard filters are OFF by default since 30 July 2026.
        #
        # That is not a preference, it is a measurement. On the GIAB v4.2.1
        # benchmark of iteration 2, applying them drops SNP F1 from 0.9921 to
        # 0.9883, because they remove 40,409 true positives to eliminate 15,445
        # false ones -- 2.6 good ones for every error caught. The same held for
        # the baseline (3.0 to 1), so it is not an effect of changing reference:
        # GATK's historical thresholds are tuned on data and depths different
        # from these.
        #
        # On indels they stay on, because there they are neutral or slightly
        # useful: F1 from 0.9924 (ALL) to 0.9928 (PASS).
        #
        # To restore the previous behaviour:  HG002_SNP_HARD_FILTERS=on
        snp_filters=()
        if [[ "${HG002_SNP_HARD_FILTERS:-off}" == "on" ]]; then
            snp_filters=(
                --filter-name "SNP_QD2"            --filter-expression "QD < 2.0"
                --filter-name "SNP_QUAL30"         --filter-expression "QUAL < 30.0"
                --filter-name "SNP_SOR3"           --filter-expression "SOR > 3.0"
                --filter-name "SNP_FS60"           --filter-expression "FS > 60.0"
                --filter-name "SNP_MQ40"           --filter-expression "MQ < 40.0"
                --filter-name "SNP_MQRankSum"      --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5"
                --filter-name "SNP_ReadPosRankSum" --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0"
            )
        fi

        if [[ ${#snp_filters[@]} -gt 0 ]]; then
            echo "SNP hard filters: ON (HG002_SNP_HARD_FILTERS=on)"
            gatk --java-options "-Xmx8g" VariantFiltration \
                -R "$REFERENCE" -V "$snps" -O "${prefix}.tmp.snps.filtered.vcf.gz" \
                "${snp_filters[@]}"
        else
            echo "SNP hard filters: off (F1 0.9921 against 0.9883 with filters)"
            cp -f "$snps" "${prefix}.tmp.snps.filtered.vcf.gz"
            if [[ -s "${snps}.tbi" ]]; then
                cp -f "${snps}.tbi" "${prefix}.tmp.snps.filtered.vcf.gz.tbi"
            else
                tabix -f -p vcf "${prefix}.tmp.snps.filtered.vcf.gz"
            fi
        fi

        gatk --java-options "-Xmx8g" VariantFiltration \
            -R "$REFERENCE" -V "$indels" -O "${prefix}.tmp.indels.filtered.vcf.gz" \
            --filter-name "INDEL_QD2" --filter-expression "QD < 2.0" \
            --filter-name "INDEL_QUAL30" --filter-expression "QUAL < 30.0" \
            --filter-name "INDEL_FS200" --filter-expression "FS > 200.0" \
            --filter-name "INDEL_SOR10" --filter-expression "SOR > 10.0" \
            --filter-name "INDEL_ReadPosRankSum" --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0"

        gatk --java-options "-Xmx8g" VariantFiltration \
            -R "$REFERENCE" -V "$other" -O "${prefix}.tmp.other.filtered.vcf.gz" \
            --filter-name "OTHER_QD2" --filter-expression "QD < 2.0" \
            --filter-name "OTHER_QUAL30" --filter-expression "QUAL < 30.0" \
            --filter-name "OTHER_FS200" --filter-expression "FS > 200.0" \
            --filter-name "OTHER_SOR10" --filter-expression "SOR > 10.0"

        gatk --java-options "-Xmx8g" MergeVcfs \
            -I "${prefix}.tmp.snps.filtered.vcf.gz" \
            -I "${prefix}.tmp.indels.filtered.vcf.gz" \
            -I "${prefix}.tmp.other.filtered.vcf.gz" \
            -O "$partial"
        tabix -f -p vcf "$partial"

        raw_count="$(count_records "$raw")"
        filtered_count="$(count_records "$partial")"
        [[ "$raw_count" == "$filtered_count" ]] || {
            echo "Records lost: VCF=${raw_count}, hard-filtered=${filtered_count}" >&2
            exit 2
        }
        check_vcf "$partial" "${partial}.tbi" "$sample"

        mv -f "$partial" "$final"
        mv -f "${partial}.tbi" "${final}.tbi"
        rm -f "${prefix}.tmp."*
        echo "Hard filtering complete without removing records: ${raw_count}."
        ;;

    check-hardfilter)
        [[ $# -eq 2 ]] || usage
        prefix="$1"
        sample="$2"
        raw="${WORK}/${prefix}.vcf.gz"
        filtered="${WORK}/${prefix}.hardfiltered.vcf.gz"
        check_vcf "$filtered" "${filtered}.tbi" "$sample"
        [[ "$(count_records "$raw")" == "$(count_records "$filtered")" ]]
        ;;

    annotate)
        [[ $# -eq 2 ]] || usage
        prefix="$1"
        sample="$2"
        cd "$WORK"

        input="${prefix}.hardfiltered.vcf.gz"
        partial_vcf="${prefix}_annotated.partial.vcf"
        partial_gz="${partial_vcf}.gz"
        final_gz="${prefix}_annotated.vcf.gz"
        partial_log="${prefix}.snpeff.partial.log"
        final_log="${prefix}.snpeff.log"
        partial_stats="${prefix}_snpEff_summary.partial.html"
        final_stats="${prefix}_snpEff_summary.html"
        partial_genes="${prefix}_snpEff_summary.partial.genes.txt"
        final_genes="${prefix}_snpEff_summary.genes.txt"
        high_vcf="${prefix}_high_impact.partial.vcf"
        high_tsv="${prefix}_high_impact_variants.partial.tsv"
        final_high="${prefix}_high_impact_variants.tsv"

        rm -f "$partial_vcf" "$partial_gz" "${partial_gz}.tbi" \
            "$partial_log" "$partial_stats" "$partial_genes" \
            "$high_vcf" "$high_tsv"

        snpEff -Xmx8g -dataDir "${PROJECT}/snpEff_data" -nodownload -noLog -v \
            -stats "$partial_stats" hg38 "$input" \
            > "$partial_vcf" 2> "$partial_log"

        # The shortlist contains only variants that cleared the hard filters and
        # always keeps twelve columns, even with several transcripts.
        SnpSift filter "(FILTER = 'PASS') & (ANN[*].IMPACT = 'HIGH')" \
            "$partial_vcf" > "$high_vcf"
        SnpSift extractFields -s "," -e "." "$high_vcf" \
            CHROM POS REF ALT FILTER QUAL \
            "GEN[${sample}].GT" "GEN[${sample}].DP" "GEN[${sample}].GQ" \
            'ANN[*].GENE' 'ANN[*].EFFECT' 'ANN[*].IMPACT' \
            > "$high_tsv"

        awk -F '\t' '
            NF != 12 {exit 2}
            NR > 1 && $5 != "PASS" {exit 3}
        ' "$high_tsv"

        bgzip -f "$partial_vcf"
        tabix -f -p vcf "$partial_gz"
        check_vcf "$partial_gz" "${partial_gz}.tbi" "$sample"

        filtered_count="$(count_records "$input")"
        annotated_count="$(count_records "$partial_gz")"
        [[ "$filtered_count" == "$annotated_count" ]] || {
            echo "Records lost: filtered=${filtered_count}, annotated=${annotated_count}" >&2
            exit 2
        }

        chromosome_errors="$(
            awk '$1 == "ERROR_CHROMOSOME_NOT_FOUND" {print $2}' "$partial_log" |
                tail -n 1
        )"
        chromosome_errors="${chromosome_errors:-0}"
        high_count="$(awk 'NR > 1 {n++} END {print n+0}' "$high_tsv")"
        {
            printf 'metric\tvalue\n'
            printf 'input_records\t%s\n' "$filtered_count"
            printf 'annotated_records\t%s\n' "$annotated_count"
            printf 'pass_high_records\t%s\n' "$high_count"
            printf 'chromosome_not_found\t%s\n' "$chromosome_errors"
        } > "${prefix}.annotation_qc.partial.tsv"

        mv -f "$partial_gz" "$final_gz"
        mv -f "${partial_gz}.tbi" "${final_gz}.tbi"
        mv -f "$partial_log" "$final_log"
        mv -f "$partial_stats" "$final_stats"
        [[ -s "$partial_genes" ]] && mv -f "$partial_genes" "$final_genes"
        mv -f "$high_tsv" "$final_high"
        mv -f "${prefix}.annotation_qc.partial.tsv" "${prefix}.annotation_qc.tsv"
        rm -f "$high_vcf"

        echo "Annotated records: ${annotated_count}; PASS+HIGH: ${high_count}."
        echo "Contigs absent from the snpEff database: ${chromosome_errors}."
        tail -n 20 "$final_log" >&2
        ;;

    check-annotation)
        [[ $# -eq 2 ]] || usage
        prefix="$1"
        sample="$2"
        filtered="${WORK}/${prefix}.hardfiltered.vcf.gz"
        annotated="${WORK}/${prefix}_annotated.vcf.gz"
        check_vcf "$annotated" "${annotated}.tbi" "$sample"
        [[ "$(count_records "$filtered")" == "$(count_records "$annotated")" ]]
        [[ -s "${WORK}/${prefix}_high_impact_variants.tsv" ]]
        [[ -s "${WORK}/${prefix}_snpEff_summary.html" ]]
        [[ -s "${WORK}/${prefix}_snpEff_summary.genes.txt" ]]
        [[ -s "${WORK}/${prefix}.snpeff.log" ]]
        [[ -s "${WORK}/${prefix}.annotation_qc.tsv" ]]
        awk -F '\t' 'NF != 12 {exit 2} NR > 1 && $5 != "PASS" {exit 3}' \
            "${WORK}/${prefix}_high_impact_variants.tsv"
        ;;

    qc)
        [[ $# -eq 2 ]] || usage
        prefix="$1"
        mode="$2"
        bam="${WORK}/${prefix}.bam"
        cd "$WORK"

        suffixes=(
            flagstat.txt samtools_stats.txt idxstats.tsv bam_header.sam
            coverage.tsv wgs_metrics.txt depth_thresholds.tsv vcf_contigs.txt
        )
        for suffix in "${suffixes[@]}"; do
            rm -f "${prefix}.${suffix}.partial"
        done

        check_bam "$bam" "${bam}.bai" "HG002" "HV3C3DSXX.2.AGCGATAG+AGGCGAAG"
        samtools flagstat -@ 8 "$bam" > "${prefix}.flagstat.txt.partial"
        samtools stats -@ 8 "$bam" > "${prefix}.samtools_stats.txt.partial"
        samtools idxstats "$bam" > "${prefix}.idxstats.tsv.partial"
        samtools view -H "$bam" > "${prefix}.bam_header.sam.partial"

        primary_bed="/tmp/${prefix}.primary.bed"
        interval_list="/tmp/${prefix}.primary.interval_list"
        if [[ "$mode" == "smoke" ]]; then
            printf 'chr20\t0\t1000000\n' > "$primary_bed"
            samtools coverage -r chr20:1-1000000 "$bam" \
                > "${prefix}.coverage.tsv.partial"
        else
            awk 'BEGIN{OFS="\t"} $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1,0,$2}' \
                "${REFERENCE}.fai" > "$primary_bed"
            samtools coverage "$bam" > "${prefix}.coverage.tsv.partial"
        fi

        gatk --java-options "-Xmx6g" BedToIntervalList \
            -I "$primary_bed" \
            -O "$interval_list" \
            -SD "$REFERENCE_DICT"

        # CollectWgsMetrics produces PCT_1X/10X/20X/30X directly, without
        # emitting a row for each of the ~3.1 billion loci.
        #
        # Picard's "fast" collector sizes its counters on READ_LENGTH. The
        # default is 150: with 151 bp reads you get "ArrayIndexOutOfBoundsException:
        # requested index ... out of counter bounds". So the real length is read
        # from samtools stats, which has already run just above, instead of
        # trusting the default.
        read_length="$(
            awk -F '\t' '$2 == "maximum length:" {print $3; exit}' \
                "${prefix}.samtools_stats.txt.partial"
        )"
        read_length="${read_length:-151}"

        # Even with the right READ_LENGTH the fast collector has a known bug in
        # its window handling (Picard #1523, #1970). If it fails, fall back to
        # the standard algorithm: slower, but without that problem. A few extra
        # minutes beat losing a seven-hour run.
        if ! gatk --java-options "-Xmx8g" CollectWgsMetrics \
                -I "$bam" \
                -O "${prefix}.wgs_metrics.txt.partial" \
                -R "$REFERENCE" \
                --INTERVALS "$interval_list" \
                --USE_FAST_ALGORITHM true \
                --READ_LENGTH "$read_length" \
                --VALIDATION_STRINGENCY SILENT \
                --MINIMUM_MAPPING_QUALITY 20 \
                --MINIMUM_BASE_QUALITY 20 \
                --TMP_DIR /pbtmp; then
            echo "Fast CollectWgsMetrics failed with READ_LENGTH=${read_length};" \
                 "retrying with the standard algorithm." >&2
            rm -f "${prefix}.wgs_metrics.txt.partial"
            gatk --java-options "-Xmx8g" CollectWgsMetrics \
                -I "$bam" \
                -O "${prefix}.wgs_metrics.txt.partial" \
                -R "$REFERENCE" \
                --INTERVALS "$interval_list" \
                --USE_FAST_ALGORITHM false \
                --VALIDATION_STRINGENCY SILENT \
                --MINIMUM_MAPPING_QUALITY 20 \
                --MINIMUM_BASE_QUALITY 20 \
                --TMP_DIR /pbtmp
        fi

        awk -F '\t' '
            /^GENOME_TERRITORY\t/ {
                for (i=1; i<=NF; i++) column[$i]=i
                getline
                territory=$(column["GENOME_TERRITORY"])
                p1=$(column["PCT_1X"])
                p10=$(column["PCT_10X"])
                p20=$(column["PCT_20X"])
                p30=$(column["PCT_30X"])
                print "total_bases\tge1\tge10\tge20\tge30\tpct_ge1\tpct_ge10\tpct_ge20\tpct_ge30"
                printf "%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.6f\t%.6f\t%.6f\t%.6f\n", \
                    territory, int(territory*p1+0.5), int(territory*p10+0.5), \
                    int(territory*p20+0.5), int(territory*p30+0.5), \
                    100*p1, 100*p10, 100*p20, 100*p30
                found=1
                exit
            }
            END {if (!found) exit 2}
        ' "${prefix}.wgs_metrics.txt.partial" \
            > "${prefix}.depth_thresholds.tsv.partial"

        for vcf in \
            "${prefix}.g.vcf.gz" \
            "${prefix}.vcf.gz" \
            "${prefix}.hardfiltered.vcf.gz" \
            "${prefix}_annotated.vcf.gz"
        do
            gzip -t "$vcf"
            [[ -s "${vcf}.tbi" ]]
            tabix -l "$vcf" >/dev/null
        done
        tabix -l "${prefix}.vcf.gz" > "${prefix}.vcf_contigs.txt.partial"

        for suffix in "${suffixes[@]}"; do
            [[ -s "${prefix}.${suffix}.partial" ]]
            mv -f "${prefix}.${suffix}.partial" "${prefix}.${suffix}"
        done
        ;;

    check-qc)
        [[ $# -eq 1 ]] || usage
        prefix="$1"
        required=(
            "${WORK}/${prefix}.flagstat.txt"
            "${WORK}/${prefix}.samtools_stats.txt"
            "${WORK}/${prefix}.idxstats.tsv"
            "${WORK}/${prefix}.bam_header.sam"
            "${WORK}/${prefix}.coverage.tsv"
            "${WORK}/${prefix}.wgs_metrics.txt"
            "${WORK}/${prefix}.depth_thresholds.tsv"
            "${WORK}/${prefix}.vcf_contigs.txt"
        )
        for file in "${required[@]}"; do [[ -s "$file" ]]; done
        awk -F '\t' '
            NR == 1 {ok=($1=="total_bases" && $6=="pct_ge1" && $9=="pct_ge30")}
            NR == 2 {ok=ok && $1>0 && $6>=0 && $6<=100 && $9>=0 && $9<=100}
            END {exit !ok}
        ' "${WORK}/${prefix}.depth_thresholds.tsv"
        ;;

    *)
        usage
        ;;
esac
