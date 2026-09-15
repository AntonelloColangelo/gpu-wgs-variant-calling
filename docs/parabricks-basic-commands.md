# Simple Parabricks commands

Run commands from the repository root in a **Bash** terminal on Linux or Ubuntu
WSL2. On WSL2, keep the working directory in the Linux filesystem, for example
`~/gpu-wgs-variant-calling`, rather than `/mnt/c/`.

## 1. Prerequisites

- NVIDIA GPU, compatible driver and Docker configured for GPU access.
- Sufficient RAM and disk space for whole-genome alignment. The historical
  experiment used a 24 GB RTX 3090 with memory-saving options and 64 GB swap;
  those measurements do not guarantee a run on another machine.
- Parabricks image pinned to the version used in the project:

```bash
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
```

The following reference-preparation commands also require `curl`, `gzip`,
`bwa`, `samtools`, `gatk` and `tabix` in the Bash environment. They are setup
tools, not dependencies of PowerShell or a local custom Docker image.

## 2. Prepare inputs

The [data sources](data-sources.md) describe the original HG002 files and
reference. Place the two HG002 FASTQ files in the repository root, using the
names shown in `run_parabricks.sh`. Verify downloads using their source checksums
when available, and check compressed-file integrity:

```bash
gzip -t HG002.novaseq.pcr-free.40x.R1.fastq.gz
gzip -t HG002.novaseq.pcr-free.40x.R2.fastq.gz
```

For a new setup, download the no-ALT reference and prepare its indexes. Run once:

```bash
mkdir -p ref
curl -fL --retry 3 -o ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gzip -t ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gunzip ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
bwa index ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
samtools faidx ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
gatk CreateSequenceDictionary -R ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna -O ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict
curl -fL --retry 3 -o ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
gzip -t ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz
tabix -p vcf ref/Homo_sapiens_assembly38.dbsnp138.vcf.gz
```

Building BWA indexes takes time. Pre-built NCBI indexes are another option,
listed in the data sources. The historical `scripts/fetch_reference_noalt.sh`
was written for the original environment and additionally requires its custom
Docker image and the old reference `.fai`; it is not the fresh-setup route above.

## 3. Read the four commands

[run_parabricks.sh](../run_parabricks.sh) contains the full commands without
helper functions, checkpoints or scheduling:

| Command | Input | Output |
|---|---|---|
| `pbrun fq2bam` | Paired FASTQ and reference | Sorted BAM with duplicates marked |
| `pbrun bqsr` | BAM, reference and known sites | Recalibration table |
| `pbrun haplotypecaller` | BAM, reference and recalibration table | gVCF |
| `pbrun genotypegvcf` | gVCF and reference | Variant VCF |

`--in-fq` receives R1, R2 and the read group together. The default read group
belongs to HG002; when changing samples, update `SAMPLE`, both FASTQ paths and
all read-group metadata, including `SM`. For multiple lanes, adapt `--in-fq`
to the actual library and lane structure rather than reusing the HG002 group.
All input paths and `OUT` are **relative to the repository root**.

The mount `$PWD:/input:ro` exposes inputs read-only inside Docker. The output
folder is mounted separately at `/output`; all `pbrun` paths use those container
paths. BQSR runs separately to release alignment memory before recalibration.
HaplotypeCaller uses the recalibration table while reading the original BAM.

```bash
bash run_parabricks.sh
```

The script runs in the foreground and creates a fresh output directory. If it
fails, inspect the error and the outputs before deciding which command to rerun.
It has no automatic resume. Keep the terminal open while it runs. A successful
exit means the tools completed, not that biological accuracy was validated.

## 4. Inspect the outputs

With `samtools` and `bcftools` available locally:

```bash
samtools quickcheck -v output/HG002_haplotypecaller/HG002.bam
samtools flagstat output/HG002_haplotypecaller/HG002.bam
bcftools stats output/HG002_haplotypecaller/HG002.vcf.gz
```

These are technical checks; the GIAB comparison below evaluates accuracy within
its defined benchmark regions.

## 5. Run the accuracy benchmark separately

`benchmark_giab.sh` needs Docker, Bash, GNU coreutils, `curl` and Python 3.
It downloads the matching GIAB truth set and runs hap.py in its own container.
The reference must match the one used for calling. The default sample is HG002:

```bash
bash benchmark_giab.sh
```

For another VCF, specify a path inside the repository and a new result label:

```bash
SAMPLE=HG003 QUERY_VCF=output/HG003_deepvariant/HG003.deepvariant.vcf.gz BENCH_LABEL=HG003_dv_simple bash benchmark_giab.sh
```

Only HG002 and HG003 are supported. Results go to `reports/giab/<BENCH_LABEL>/`.
An existing directory is rejected to preserve earlier results. Set a new label
to repeat a benchmark. The new core script creates an unfiltered VCF, so its
`PASS` results are not directly equivalent to the historical filtered call set.

## 6. DeepVariant on the same alignment

Edit `SAMPLE` and `BAM` in
[run_deepvariant.sh](../hg003-deepvariant-vs-gatk/run_deepvariant.sh), then run:

```bash
bash hg003-deepvariant-vs-gatk/run_deepvariant.sh
```

The default is HG003. Use the sorted, duplicate-marked BAM with original base
qualities and its index. This script makes one `pbrun deepvariant` call and
writes to a separate output folder. It does not require `config/sample.env`.

## 7. Historical post-processing

`scripts/postprocess.sh` retains the filtering, snpEff annotation and QC code
used for the original results. It expects `/project` and `/work` mounts, GATK,
samtools, snpEff and other tools from the original `bioinfo-codeserver:latest`
image, plus a prepared snpEff database. Its image recipe is not published here.
`scripts/generate_report.py` expects the original collection of QC files and
logs. They are retained as analysis evidence, not automatically run by the
simple script. Neither the simplified scripts nor these historical dependencies
have been executed as part of this repository cleanup.

## Documentation

- [NVIDIA Parabricks 4.7.0 tutorials](https://docs.nvidia.com/clara/parabricks/4.7.0/tutorials.html)
- [NVIDIA whole-genome calling example](https://docs.nvidia.com/clara/parabricks/latest/tutorials/how-tos/wholegenomegermlinesmallvariants.html)
