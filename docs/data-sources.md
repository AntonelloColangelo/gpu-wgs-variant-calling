# Data sources

The genomic inputs are public but too large for GitHub. Download them into the
paths shown below before running `./run_hg002.sh`.

## GRCh38 no-ALT + hs38d1 reference

The final run used NCBI's no-ALT analysis set. It retains the hs38d1 decoys but
removes ALT/HLA scaffolds, avoiding the ALT-without-`.alt` mismatch described in
the README.

Base directory:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
```

Download these files into `ref/`:

```text
GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai
GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz
md5checksums.txt
```

Decompress the FASTA, unpack its five BWA index files beside it, and retain the
published `.fai`. The filename expected by the script is:

```text
ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
```

## HG002 FASTQ

NovaSeq PCR-free paired-end reads distributed as HG002 40× by Google Brain
Genomics:

```text
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG002.novaseq.pcr-free.40x.R1.fastq.gz
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG002.novaseq.pcr-free.40x.R2.fastq.gz
```

Expected sizes are 34,156,056,301 bytes (R1) and 35,449,847,992 bytes (R2),
containing 474,384,500 read pairs. Store both files in the repository root.

## dbSNP known sites

The BQSR step uses Broad's GRCh38 dbSNP 138 VCF and its index:

```text
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
```

Store them in `ref/`.

## GIAB benchmark truth set

The published accuracy figures use HG002 GRCh38 v4.2.1 and are restricted to
the high-confidence BED:

```text
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

These files are only needed to reproduce the hap.py benchmark, not to run the
variant-calling pipeline.
