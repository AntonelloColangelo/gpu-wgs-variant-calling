# Whole-genome variant calling with NVIDIA Parabricks

An independent bioinformatics project using a single RTX 3090 to analyse
human whole-genome sequencing data. Experiments use the public GIAB reference
samples **HG002** and **HG003**.

The main code is now a short Bash script with four explicit Parabricks commands:

```text
paired FASTQ → fq2bam → BAM → bqsr → HaplotypeCaller → gVCF → GenotypeGVCF → VCF
```

**[Read the script](run_parabricks.sh)** ·
**[Setup and commands](docs/parabricks-basic-commands.md)** ·
**[HG003 caller comparison](hg003-deepvariant-vs-gatk/README.md)**

## Run the analysis

Use a Bash terminal in Linux or Ubuntu WSL2, with Docker and NVIDIA GPU access.
Work in a Linux filesystem (under the WSL home directory on Windows), and run
commands from the repository root. The script uses Parabricks **4.7.0-1**, the
version used in the documented experiments.

1. Prepare the reference, its indexes, paired FASTQ files and dbSNP known sites
   using the [setup guide](docs/parabricks-basic-commands.md).
2. Edit `SAMPLE`, `R1`, `R2` and `READ_GROUP` at the top of
   [run_parabricks.sh](run_parabricks.sh). The defaults describe HG002.
3. Run:

```bash
bash run_parabricks.sh
```

The main outputs are `output/HG002_haplotypecaller/HG002.bam`,
`HG002.recal.txt`, `HG002.g.vcf.gz` and `HG002.vcf.gz` in the same directory.
The final VCF contains variant calls before downstream filtering and annotation.

The commands run in the foreground, stop on error and require a fresh output
directory. There is no scheduler, checkpoint recovery or background monitoring.
For another run, choose a new directory, for example:

```bash
OUT=output/HG002_run2 bash run_parabricks.sh
```

For the separate DeepVariant experiment on an existing HG003 BAM:

```bash
bash hg003-deepvariant-vs-gatk/run_deepvariant.sh
```

Set its `SAMPLE` and `BAM` variables to match your actual alignment. Both callers
use the same alignment; HaplotypeCaller applies the BQSR table during calling,
while DeepVariant reads the BAM with its original qualities.

## Published experiments

These measurements come from the original workflow. They have **not** been
remeasured with the simplified scripts. The original code remains available
in [Git history](https://github.com/AntonelloColangelo/gpu-wgs-variant-calling/tree/f222fdf).

### Accuracy against Genome in a Bottle

`hap.py` with the `vcfeval` engine, GIAB v4.2.1, GRCh38, chr1–22 within the
high-confidence BED. Values below use the `ALL` rows.

| Sample | Caller | SNP F1 | INDEL F1 | Benchmark results |
|---|---|---:|---:|---|
| HG002 | Parabricks HaplotypeCaller | 0.9921 | 0.9924 | [hap.py summary (CSV)](reports/giab/HG002_NovaSeq_40x_53007e55/happy.summary.csv) |
| HG003 | Parabricks HaplotypeCaller | 0.9917 | 0.9913 | [hap.py summary (CSV)](hg003-deepvariant-vs-gatk/giab/haplotypecaller/happy.summary.csv) |
| HG003 | Parabricks DeepVariant | 0.9962 | 0.9959 | [hap.py summary (CSV)](hg003-deepvariant-vs-gatk/giab/deepvariant/happy.summary.csv) |

The linked CSV files are the original `hap.py` summaries, including true-positive,
false-positive and false-negative counts, precision, recall and F1 scores.
[benchmark_giab.sh](benchmark_giab.sh) runs the accuracy comparison separately;
see the setup guide for selecting the sample and VCF.

### Hardware and performance

The historical runs used an RTX 3090 (24 GB) connected through OCuLink to a
Lenovo Legion 5 with a Ryzen 5 6600H, 32 GB RAM and Ubuntu WSL2.
HG002 took **10 h 31 min** through the original full workflow, including
post-processing. This is not the runtime of the new four-command script.

A separate CPU/GPU experiment measured alignment, sorting and duplicate marking
on 10 million read pairs, and variant calling on `chr20:1-20000000` at full
coverage. GPU speedups on those workloads were approximately **7.3×** and
**4.7×**, respectively. See the [execution times (TSV)](reports/cpu_vs_gpu/benchmark.tsv)
and [benchmark code](scripts/benchmark_cpu_vs_gpu.sh); full-genome CPU time was
not measured. That historical benchmark script requires its original Docker
volumes and local tools image.

## Analysis notes and scope

- The reference is GRCh38 **no-ALT + hs38d1 decoys**. The
  [historical analysis](docs/historical-analysis.md) documents the original
  reference problem, the subsequent improvement and the resource measurements.
- BQSR uses dbSNP138, not a complete modern known-sites bundle.
- Filtering, snpEff annotation, QC and HTML report code from the original analysis
  is retained in `scripts/postprocess.sh` and `scripts/generate_report.py`.
  These are separate, historical tools with additional dependencies; the simple
  script stops at the variant VCF. See the setup guide before using them.
- Accuracy figures apply to the benchmark regions, not the entire genome or
  clinical interpretation. Functional impact annotations are not pathogenicity
  classifications.
- FASTQ, BAM, VCF, references and databases are excluded from Git. Only code,
  documentation and small result files are published.

## Repository map

| File or directory | Purpose |
|---|---|
| `run_parabricks.sh` | Four linear Parabricks commands, FASTQ to VCF |
| `hg003-deepvariant-vs-gatk/run_deepvariant.sh` | One DeepVariant command on an existing BAM |
| `benchmark_giab.sh` | Accuracy comparison against HG002 or HG003 GIAB truth |
| `scripts/` | Reference preparation and historical post-processing/benchmark tools |
| `reports/` | Saved HG002 accuracy, timing and smoke-test results |
| `hg003-deepvariant-vs-gatk/` | HG003 comparison, report and saved accuracy results |
| `docs/` | Setup, data sources and historical experiment notes |

## References

- [NVIDIA Parabricks documentation](https://docs.nvidia.com/clara/parabricks/4.7.0/tutorials.html)
- [Input data sources](docs/data-sources.md)
- [Historical HG002 analysis](docs/historical-analysis.md)
