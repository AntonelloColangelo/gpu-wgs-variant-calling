# 1) FQ2BAM

Teaching notes on the two main steps. The commands below run against the reduced
FASTQ files in `smoke/` (one million read pairs): they exist to show the
mechanics, not to produce the final result. The real analysis is launched with
`.\Start-HG002.ps1`, which adds BQSR, hard filtering, annotation and integrity
checks.

Run these from Ubuntu WSL2. Change into the project directory first, so that
`$(pwd)` expands to `/mnt/c/HG002 Parabricks Experiment`.

In the `--volume` flags, the path to the left of the colon is the host path
(`/mnt/c/...`); the one on the right exists only inside the container. From that
point on `pbrun` must be given `/workdir` and `/outputdir`, never `/mnt/c/...`.
The project is mounted read-only (`:ro`) so the reference and the FASTQ files
cannot be touched by accident.

```bash
cd "/mnt/c/HG002 Parabricks Experiment"
mkdir -p output/tmp

docker run \
    --gpus all \
    --rm \
    --volume "$(pwd)":/workdir:ro \
    --volume "$(pwd)/output":/outputdir \
  nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1 \
  pbrun fq2bam \
    --ref /workdir/ref/Homo_sapiens_assembly38.fasta \
    --in-fq /workdir/smoke/HG002_NovaSeq_smoke_R1.fastq.gz /workdir/smoke/HG002_NovaSeq_smoke_R2.fastq.gz '@RG\tID:HV3C3DSXX.2\tPL:ILLUMINA\tPM:NovaSeq6000\tLB:HG002_PCR_FREE\tPU:HV3C3DSXX.2.AGCGATAG+AGGCGAAG\tSM:HG002' \
    --out-bam /outputdir/fq2bam_output.bam \
    --tmp-dir /outputdir/tmp \
    --low-memory
```

Both FASTQ files go on the same `--in-fq` flag: given only one, Parabricks
treats the run as single-end. The third element is the read group — without it
the BAM comes out with a generic `@RG` and the sample is not named `HG002`.

If you get an out-of-memory error make sure your computer has enough RAM, and
that large amounts of memory aren't being used by other programs.


# 2) HAPLOTYPECALLER

```bash
cd "/mnt/c/HG002 Parabricks Experiment"

docker run \
    --gpus all \
    --rm \
    --volume "$(pwd)":/workdir:ro \
    --volume "$(pwd)/output":/outputdir \
  nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1 \
  pbrun haplotypecaller \
    --ref /workdir/ref/Homo_sapiens_assembly38.fasta \
    --in-bam /outputdir/fq2bam_output.bam \
    --out-variants /outputdir/variants.vcf \
    -L chr20:1-1000000 \
    --tmp-dir /outputdir/tmp \
    --htvc-low-memory
```

In HaplotypeCaller the flag is called `--htvc-low-memory`; `--low-memory` exists
only for fq2bam. `-L` restricts calling to a single window: without it,
HaplotypeCaller scans the whole genome even when there are very few reads. Note
that `-L` does not exist for fq2bam — alignment is always genome-wide, so
restricting regions saves no time on the expensive step.

The outputs are .bam, .bam.bai, .txt and .vcf files.
The VCF with the variants should show lines like this:

```
chr1    16378   .   T   C   45.28   .   AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=23.55;QD=22.64;SOR=2.303    GT:AD:DP:GQ:PL  1/1:0,2:2:6:57,6,0
```
