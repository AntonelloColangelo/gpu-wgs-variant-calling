# DeepVariant vs. GATK HaplotypeCaller on HG003 — a clean rematch

**Sample**: HG003 / NA24149, GIAB Ashkenazim trio (father), Illumina NovaSeq PCR-free 35x
**Same BAM, same truth set, same BED — the only variable is the variant caller.**
**Full interactive report (Italian)**: [`report.html`](report.html)

The historical experiment used the same germline pipeline on a second GIAB
sample. The simplified command is now published in [`run_deepvariant.sh`](run_deepvariant.sh)
calls DeepVariant on the same BAM HaplotypeCaller already produced, so the
two callers are compared on identical input. Raw hap.py evidence for both
runs is under [`giab/`](giab).

## Why HG003, and not HG002 again

This project already evaluated DeepVariant once, for the main HG002 run, and
dropped the comparison: Google's public DeepVariant WGS model is trained on
GIAB samples HG001–HG007 *except HG003*. Benchmarking it against HG002 would
have graded it on data it had already seen. HG003 is the one genome in the
trio it never trained on — the only fair ground for this test.

## Accuracy against GIAB v4.2.1 (hap.py / vcfeval, chr1–22 high-confidence BED)

| Metric | HaplotypeCaller | DeepVariant | Δ |
|---|---:|---:|---:|
| SNP recall | 0.9922 | 0.9937 | +0.0015 |
| SNP precision | 0.9912 | 0.9986 | +0.0074 |
| SNP F1 | 0.9917 | 0.9962 | +0.0045 |
| INDEL recall | 0.9906 | 0.9942 | +0.0036 |
| INDEL precision | 0.9919 | 0.9976 | +0.0057 |
| INDEL F1 | 0.9913 | 0.9959 | +0.0046 |

Both classes improve on both axes at once — usually precision gains cost
sensitivity, not here. Full numbers: [`giab/haplotypecaller/giab_benchmark.json`](giab/haplotypecaller/giab_benchmark.json),
[`giab/deepvariant/giab_benchmark.json`](giab/deepvariant/giab_benchmark.json).

![F1 for SNP and INDEL, HaplotypeCaller vs DeepVariant](figures/fig1_f1.png)

| Errors (absolute) | HaplotypeCaller | DeepVariant | Change |
|---|---:|---:|---:|
| SNP false positives | 29,393 | 4,567 | −84.5% |
| SNP false negatives | 25,832 | 20,985 | −18.8% |
| INDEL false positives | 4,243 | 1,280 | −69.8% |
| INDEL false negatives | 4,744 | 2,932 | −38.2% |

![False positive counts for SNP and INDEL, HaplotypeCaller vs DeepVariant](figures/fig2_false_positives.png)

DeepVariant reports ~230,000 fewer `PASS` variants overall, which looks like
a red flag until you check it against ground truth: inside the judged
regions it removes ~25,000 wrong SNPs and adds ~4,800 right ones. Two
independent, non-F1 indicators agree with the direction: GIAB's known
het/hom ratio is 1.535 (HaplotypeCaller 1.670, DeepVariant 1.447), and the
expected Ti/Tv is 2.103 (HaplotypeCaller 1.962, DeepVariant 1.979).

## Runtime

Variant calling: **86 min (DeepVariant) vs. 4h 24m (HaplotypeCaller)** on the
same RTX 3090 — about **3x faster**, not the predicted outcome since
DeepVariant does more per-site work (pileup images + CNN inference). The gain
comes from GPU parallelism that HaplotypeCaller barely uses. DeepVariant also
needs no BQSR and no hard filtering, removing two pipeline stages worth ~18
more minutes.

![Variant calling time, HaplotypeCaller vs DeepVariant](figures/fig3_runtime.png)

## The clinical stress test

ClinVar annotation surfaced 7 filter-surviving Pathogenic/Likely pathogenic
calls. Checking each against DeepVariant's independent call:

- **3 are artifacts** — DeepVariant marks them `RefCall`.
- **2 are real but too common to be pathogenic** — CDKN2B (79% population
  frequency), XG (37%).
- **2 are genuine heterozygous carrier states**, consistent with the donor's
  Ashkenazi origin.

One artifact, **PRSS1** (chr7:142750561 C>T), was reported as a confirmed
pathogenic variant by an online source with matching rsID and protein
change. `MQRankSum = −6.02` is the text signature of systematic mismapping
in the ALT-supporting reads — stronger than the same signal on SLC9B1, which
that source correctly flagged as false. The region also overlaps the TRB
locus, which is somatically rearranged in blood-derived DNA. None of this
project's active hard filters check MQ or MQRankSum on indels; DeepVariant
rejects the call independently, without needing to.

![Outcome of the 7 filter-surviving Pathogenic/Likely-pathogenic calls](figures/fig4_clinical_variants.png)

## Takeaways

- DeepVariant replaces HaplotypeCaller with no tradeoff here: better recall
  *and* precision on both variant classes, 3x faster, same hardware.
- Hard filtering and VQSR become moot: DeepVariant's unfiltered precision
  (0.9986) beats HaplotypeCaller's best filtered result. Its `ALL` and
  `PASS` benchmark rows are identical.
- The benchmark still only judges chr1–22 inside the high-confidence BED:
  `Frac_NA` shows 13% of called SNPs and 44% of called indels aren't judged
  at all, and chrX/chrY/chrM have no truth set here.
- No population-frequency filter yet: without gnomAD, common variants like
  CDKN2B and XG still enter the clinical table.

---
*Truth set NIST/GIAB v4.2.1 for HG003 (NA24149), GRCh38, chr1–22
high-confidence regions. Reference GCA_000001405.15, no ALT contigs, hs38d1
decoy. ClinVar release 2026-07-28. HG003 is a public reference sample;
reported variants are not a clinical interpretation.*

## Running the simplified commands

Run from the repository root. Prepare HG003 FASTQ and its actual read group,
then edit those inputs in `run_parabricks.sh` and run it to produce the BAM and
HaplotypeCaller VCF. The defaults in that script describe HG002, not HG003.
Run DeepVariant on the HG003 BAM with:

```bash
bash hg003-deepvariant-vs-gatk/run_deepvariant.sh
```

Compare both VCFs using the matching HG003 truth set:

```bash
SAMPLE=HG003 BENCH_LABEL=HG003_hc_simple bash benchmark_giab.sh
SAMPLE=HG003 QUERY_VCF=output/HG003_deepvariant/HG003.deepvariant.vcf.gz BENCH_LABEL=HG003_dv_simple bash benchmark_giab.sh
```

These commands use fresh output directories. The simplified scripts have not
been rerun on a whole genome. The tables and HTML report above are historical
results, and references to the old automation in that report describe the
original experiment. Individual clinical assertions in the historical report
are exploratory and have not been independently validated by this code cleanup;
caller agreement or disagreement alone does not establish biological truth.
