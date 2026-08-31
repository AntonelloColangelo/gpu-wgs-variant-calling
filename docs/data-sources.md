# Data sources

Direct download addresses for everything the pipeline expects to find locally.
None of these files are versioned in this repository: they amount to tens of
gigabytes and all come from public archives.

## GRCh38 reference — no-ALT analysis set (current default)

Used from iteration 2 onwards. `scripts/fetch_reference_noalt.sh` downloads,
verifies and prepares it; the addresses are recorded here for provenance.

```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/md5checksums.txt
```

The BWA indexes are published pre-built, so `bwa index` is never run locally:
it would cost about an hour of CPU for an identical result. The `.dict` is
generated with `samtools dict`. Integrity is checked three ways: the MD5
declared by NCBI, a `gzip -t` on the archives, and a regenerated `.fai`
compared byte for byte against the published one.

**Why this variant and not another.** NCBI publishes four analysis sets. This
one differs from the iteration-1 reference in exactly one respect: the ALT
contigs (and the HLA contigs that depend on them) are gone, while the hs38d1
decoys are kept. The decoys matter — they absorb junk reads that would
otherwise stick to real chromosomes and produce false positives. Dropping them
together with the ALT contigs would move two variables at once.

Contig names stay UCSC-style (`chr1`, `chr2`, …) and the primary chromosomes
keep identical names and lengths, which is why dbSNP, snpEff and the GIAB BED
remain usable without modification. `scripts/fetch_reference_noalt.sh` verifies
this and refuses to finish if it is not true.

## GRCh38 reference — Broad assembly38 (iteration 1, kept as baseline)

Public GATK hg38 bundle on Google Cloud, the same one used by GATK workflows.

```
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
```

Pre-built BWA indexes, so they do not have to be regenerated locally:

```
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt
```

> **Warning.** The `.alt` file above is listed here but was never downloaded,
> and the local BWA indexes are not the `.64.*` variants either. That omission
> is the root cause of the MHC sensitivity loss documented in the README: the
> reference keeps its 261 ALT contigs while BWA has no way to know they are
> alternates.
>
> **Resolved in iteration 2** by switching to the no-ALT analysis set above,
> not by adding the `.alt` file. The preflight now refuses to start on a
> reference that carries `_alt` contigs without an `.alt` file beside it, which
> is the check that would have caught this at minute zero instead of after the
> benchmark. To reproduce the defective baseline deliberately:
>
> ```bash
> HG002_REF_NAME=Homo_sapiens_assembly38.fasta \
> HG002_ALLOW_ALT_WITHOUT_ALT_FILE=1 \
>   bash run_parabricks_hg002.sh --preflight
> ```

```
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa
```

All of these belong in `ref/`.

## HG002 FASTQ

Whole-genome NovaSeq PCR-free, paired-end, nominal 40× depth.

Measured properties of the files in use:

| Quantity | Value |
|---|---|
| Read pairs | 474,384,500 |
| Read length | 151 bp, uniform |
| Total bases | ~143.3 Gbp |
| Raw coverage over 3.1 Gbp | ~46× |
| Coverage after duplicates (12.08%) | ~41× |
| File sizes | 34.2 GB (R1) and 35.4 GB (R2) |

Source, from the public Google Brain Genomics bucket:

```
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG002.novaseq.pcr-free.40x.R1.fastq.gz
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG002.novaseq.pcr-free.40x.R2.fastq.gz
```

Provenance verified by direct comparison: the remote objects report
34,156,056,301 and 35,449,847,992 bytes, byte-for-byte identical to the local
files. "40×" is the name under which the publisher distributes the dataset; the
raw coverage computed above is higher because the nominal figure refers to the
usable coverage expected after mapping and duplicate removal.

> **Do not confuse this with the UCSC variant.** The Giraffe directory of the
> UCSC project hosts an HG002 NovaSeq PCR-free set with a deceptively similar
> name — 35×, different capitalisation (`HG002.NovaSeq.pcr-free.35x.*`), and
> sizes of 29,911,262,422 and 31,043,207,981 bytes. It is not the source used
> here.
>
> ```
> https://cgl.gi.ucsc.edu/data/giraffe/mapping/reads/real/HG002/HG002.NovaSeq.pcr-free.35x.R1.fastq.gz
> https://cgl.gi.ucsc.edu/data/giraffe/mapping/reads/real/HG002/HG002.NovaSeq.pcr-free.35x.R2.fastq.gz
> ```

## GIAB truth set for benchmarking

GIAB/NIST HG002 GRCh38 v4.2.1, the ground truth used to measure precision and
recall of the variant calling.

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

The BED delimits the high-confidence regions: the comparison must be restricted
to those, otherwise regions where not even GIAB knows the answer get counted as
errors.

## dbSNP for BQSR

The project uses `Homo_sapiens_assembly38.dbsnp138.vcf.gz`, available in the
same Broad bundle. As stated in the README limitations, this is not the complete
known-sites set from the Best Practices.
