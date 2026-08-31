# Validation of the Bash pipeline

Date: 25 July 2026  
Sample: HG002  
Mode tested: smoke test on 1,000,000 NovaSeq read pairs

## Outcome

`PASS`: the pipeline completed every step automatically with a final exit
code of 0.

| Step | Duration | Result |
|---|---:|---|
| fq2bam | 106 s | PASS |
| BAM check | 1 s | PASS |
| separate BQSR | 82 s | PASS |
| BQSR check | 1 s | PASS |
| HaplotypeCaller gVCF | 50 s | PASS |
| gVCF check | 1 s | PASS |
| GenotypeGVCFs | 3 s | PASS |
| VCF check | 1 s | PASS |
| hard filtering | 21 s | PASS |
| snpEff annotation | 55 s | PASS |
| WGS QC | 42 s | PASS |
| export | 6 s | PASS |
| report | 1 s | PASS |

## Checks passed

- BAM readable, valid EOF block, usable index;
- coordinate sorted;
- read group `SM:HG002` and `PU:HV3C3DSXX.2.AGCGATAG+AGGCGAAG`;
- BQSR report with a valid GATK structure;
- gVCF, VCF, hard-filtered VCF and annotated VCF all readable and indexed;
- 12 records in the raw VCF, 12 after hard filtering and 12 after snpEff;
- no record dropped by the filters;
- HIGH-impact table restricted to `PASS` records, with a constant 12 columns;
- WGS metrics produced with `GATK CollectWgsMetrics`;
- HTML report and JSON summary generated only once every mandatory output was
  present.

## Resume and corruption

A second launch recognised the checkpoints and skipped every bioinformatic
computation that had already been validated.

On deliberately truncated throwaway copies:

- the BAM was rejected with exit code 16;
- the compressed VCF was rejected with exit code 1.

So the resume logic does not accept a file merely because it exists.

## 40× preflight

The final preflight of full mode returned exit code 0:

- complete FASTQ files present;
- historical full check reused: 474,384,500 read pairs, `PASS`;
- reference and all indexes present;
- dbSNP and its index present;
- snpEff database present;
- RTX 3090 and Docker images available;
- roughly 908 GiB free on the Linux Docker disk;
- roughly 869 GiB free on the Windows disk.

The preflight did not start fq2bam, nor any analysis of the 40× FASTQ files.
