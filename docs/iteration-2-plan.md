# Iteration 2 plan — HG002 / Parabricks

A working document. It describes every phase of the second iteration: what each
one covers, what it costs, and how you check that it worked.

Written on 29 July 2026. Reference baseline: the complete run of 27–28 July 2026
(`logs/runs/20260727_165638_full` and `logs/runs/20260728_090009_full`), whose
results are in `output/`, `reports/HG002_NovaSeq_40x_summary.json` and
`reports/giab/`.

---

## 1. Where this starts

The first iteration is **complete and valid**: nine steps out of nine `PASS` over
474,384,500 read pairs, and a benchmark against NIST HG002 v4.2.1 that produced
F1 0.9837 on indels and 0.9814 on SNPs inside the high-confidence BED.

The same benchmark, however, exposed three real defects, all documented in
`README.md`.

### Defect 1 — ALT contigs without the `.alt` file (the most serious)

The reference `Homo_sapiens_assembly38.fasta` contains **3,366 contigs**, of
which **261 end in `_alt`** and **525 are HLA contigs**. Next to the FASTA there
is **no** `Homo_sapiens_assembly38.fasta.alt`, and the BAM header carries no `AH`
tag.

Without that file BWA has no way of knowing that an ALT contig and the
corresponding primary locus describe the same piece of genome: it treats the
double alignment as ordinary multi-mapping and assigns **MAPQ 0**.
HaplotypeCaller discards everything below MAPQ 20. The variants underneath are
never called.

Measured consequences:

| Measure | Value |
|---|---|
| Reads at MAPQ 0 | 61.4 million (6.5 %) |
| True SNPs lost, whole genome | 101,963 (3.03 %) |
| True SNPs lost on chr6 | 23,093 (**10.39 %**) |
| True SNPs lost in `chr6:28,510,120-33,480,577` (MHC) | 19,876 of 20,177 (**98.51 %**) |
| Mean depth inside the MHC | 7.64× against 35× genome-wide |
| Share of chr6 false negatives due to the MHC alone | 86.1 % |

The diagnostic signal is the inverted ordering of the recalls: 97.36 % on indels
against 96.97 % on SNPs. SNPs are the easy case. If they do worse than indels the
problem is not the variant caller: it is that the data never arrives.

### Defect 2 — the SNP hard filters make F1 worse

Going from `ALL` to `PASS` on SNPs loses **41,350 true positives** in order to
remove **13,583 false ones**: three good variants thrown away for every error
caught. F1 falls from 0.9814 to 0.9770. On indels the same filters are neutral
(0.9837 → 0.9839).

The dominant filter is `SNP_MQ40`, that is, a threshold on mapping quality —
**the very quantity that defect 1 artificially depresses**. The two defects
compound: the reference lowers MAPQ, then the filter deletes whatever survived.

### Defect 3 — diploid ploidy on chrX, chrY and chrM

HaplotypeCaller was run with a single ploidy setting, diploid, across the whole
genome. HG002 is male: outside the pseudoautosomal regions chrX and chrY are
haploid. The VCF contains **10,740 PASS records on chrY, of which 7,777 are
heterozygous (72.41 %)**, which is biologically impossible. chrM is represented
diploid, and that is not a heteroplasmy analysis under any definition.

### A false alarm already ruled out — the excess of indels

The SNP/indel ratio of the call set (4.01) is lower than that of the truth set
(6.40), which looks like an excess of indels. Inside the high-confidence regions,
though, the ratio is **6.13 against 6.40**: near normal. The excess sits
**outside** those regions — 41.8 % of the indel calls fall where GIAB does not
judge, against 13.7 % of the SNPs.

Those zones are repeats, homopolymers, alt contigs and decoys: GIAB excludes them
precisely because they are hard. The 6.40 of the truth set **is not the real
biological ratio**, it is the ratio in a genome cleaned of the regions where
indels are abundant. A share of GATK indel false positives in homopolymers
remains, but it is not the main defect and does not justify changing variant
caller.

---

## 2. Two rules that hold for every phase

**One variable at a time.** Iteration 2 changes the reference and nothing else.
Filters, BQSR, snpEff, known-sites and parameters all stay identical. That way
every difference in the results can be attributed with certainty to a single
cause. The filters get re-discussed afterwards, on the new data, without
re-running anything.

**The baseline is untouchable.** The run of 27–28 July is the term of comparison.
It must not be overwritten, moved or "updated". That applies to `output/`, to
`reports/` and to the Docker volume `hg002_work_v1`.

---

## 3. The phases

### Phase 0 — Checkpoint provenance

**Why it blocks everything else.** The resume logic that saved seven hours after
the QC crash is a trap at the next iteration. Checkpoints are indexed on a
`PREFIX` and a `WORK_VOLUME` written by hand in the script. Validation checks that
a file is structurally sound, indexed, and carries the right sample and read
group — but nothing ties a checkpoint to the identity of the FASTQ files, the
reference, the parameters or the container image.

Change the reference without touching prefix and volume and the old BAM **would
pass validation**, alignment would be skipped, and the pipeline would report
`PASS` while returning exactly the results it was supposed to replace.

**What it covers.**

1. A `compute_run_key()` function that derives an eight-character hexadecimal key
   from the ingredients that determine the result: Parabricks image ID, SHA-256
   of the reference `.fai`, name and size of both FASTQ files, name and size of
   the known-sites file, read group, calling intervals. The `.fai` is the key
   choice: a few kB, and it changes if even a single contig of the reference
   changes.
2. `PREFIX`, `WORK_VOLUME` and `TMP_VOLUME` derived from that key. Changing the
   reference automatically produces a new prefix and new volumes: the old
   checkpoints are no longer found and at the same time **are not destroyed**.
3. A plain-text `run_manifest.json` in the run directory, listing every
   ingredient of the key. A hash on its own is not verifiable: you need to be
   able to read what went into it.
4. The reference becomes a variable (`REF_NAME`, overridable with
   `HG002_REF_NAME`) instead of a path written by hand in four different files.
   The preflight's list of required files is derived from it.
5. The FASTQ integrity check stays reusable: it depends only on the FASTQ files,
   not on the reference, and redoing it would cost hours for nothing.

**Files touched.** `run_parabricks_hg002.sh`, `scripts/pipeline_functions.sh`,
`scripts/postprocess.sh`, `benchmark_giab.sh`.

**Cost.** No GPU compute.

**Acceptance criterion.** With the old reference the key reproduces a stable
prefix across two consecutive invocations; changing `HG002_REF_NAME` changes the
key; `--preflight` passes in both cases.

---

### Phase 1 — Reference without ALT contigs

**What it covers.**

1. Download from NCBI of
   `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`, of its
   published `.fai` and of the **pre-built BWA indexes** (`bwa index` over 3 Gb
   would cost about an hour of CPU for an identical result).
2. MD5 verification against NCBI's `md5checksums.txt`.
3. Regeneration of the `.fai` with `samtools faidx` and comparison with the
   published one: if they match, the FASTA is intact line by line.
4. Creation of the `.dict` dictionary, which GATK and Parabricks require.
5. **A contig diff** against the current reference, with an explicit check that
   the primary chromosomes change neither name nor length — that is the condition
   under which dbSNP, snpEff and the GIAB BED stay compatible without changes.

**Why the `no_alt_plus_hs38d1` variant.** There are four of them. This one changes
a single thing with respect to today: it removes the ALT contigs and the HLA
contigs that depend on them, and **keeps the decoys**. The decoys earn their
place: they are a bin that absorbs junk reads which would otherwise stick to the
real chromosomes and create false positives. The current reference has them, and
removing them alongside the ALTs would mean moving two variables instead of one.

**Cost.** About 4 GB of download. No GPU compute.

**Acceptance criterion.** Contigs from 3,366 to 2,580, with 786 removed
(261 `_alt` + 525 HLA) and none added; lengths of chr1, chr6, chr20, chrX, chrY
and chrM unchanged.

**Script.** `scripts/fetch_reference_noalt.sh`, re-runnable and idempotent.

---

### Phase 2 — Smoke test on the new reference

**Why.** To find out in a quarter of an hour, rather than after ten, whether the
`.dict` is wrong, whether the known-sites file does not match, whether GATK
rejects the `.fna` extension, whether snpEff complains. That is exactly why the
smoke test exists.

**What it covers.** The complete nine-step pipeline over one million read pairs
(`smoke/HG002_NovaSeq_smoke_R1/R2.fastq.gz`), with calling restricted to
`chr20:1-1,000,000`. New volume and new prefix, derived from the phase 0 key.

**Cost.** About 15 minutes.

**Acceptance criterion.** Nine steps `PASS`, BAM and VCF validated, HTML report
generated, and no reuse of checkpoints belonging to the old reference.

---

### Phase 3 — Full run

**What it covers.** The pipeline over the full 474,384,500 read pairs, with the
new reference and **everything else identical**: same hard filters, same BQSR,
same flags `--bwa-options=-Y --low-memory --memory-limit 8
--bwa-normalized-queue-capacity 2 --gpuwrite`.

**Cost.** About 10 h 40 m, split as follows according to the baseline: fq2bam
2 h 01, BQSR 17 m, HaplotypeCaller 4 h 41, GenotypeGVCFs 1 m 37, hard filtering
3 m 12, snpEff 19 m, WGS QC 2 h 40, export 30 m, report 1 m 41.

**Space.** 612 GB free on `C:` and 952 GB on the Linux disk: the new BAM takes
about 86 GB and the baseline has to be preserved. There is margin, but the WSL
`.vhdx` grows in parallel and needs watching.

**Acceptance criterion.** Nine steps `PASS`, and a key in the run manifest
different from the baseline's.

---

### Phase 4 — Benchmark and verification of the predictions

**What it covers.**

1. `bash benchmark_giab.sh` on the new VCF, against the same NIST v4.2.1 truth
   set and the same high-confidence BED.
2. Recomputation of false negatives per chromosome.
3. A targeted recomputation over `chr6:28,510,120-33,480,577`, the point where
   the diagnosis is either confirmed or falls apart.
4. A count of MAPQ 0 reads in the new BAM.

**The predictions have to be written before measuring.** These are they, and they
are inferences, not measurements:

| Metric | Baseline | Expected |
|---|---|---|
| Reads at MAPQ 0, whole genome | 6.49 % | ~4.9 % — **prediction corrected on 29/07, see below** |
| Reads at MAPQ 0 **inside the MHC** | 91.45 % | 1–2 % |
| Depth inside the MHC | 7.64× | ~35×, in line with the genome |
| SNP recall inside the MHC | 1.49 % | > 95 % |
| SNPs lost on chr6 | 10.39 % | ~2 %, in line with the rest |
| SNP recall genome-wide | 96.97 % | ~99 % |
| SNP F1 (ALL) | 0.9814 | 0.993–0.995 |
| Ordering of the recalls | indel > SNP | SNP > indel |

> **Corrected prediction.** The first draft of this plan expected MAPQ 0 reads to
> fall from 6.5 % to 1–2 % **across the whole genome**. The phase 2 smoke test
> shows 4.86 %, and that prediction was wrong. The reason: ALT contigs cover a
> small fraction of GRCh38, so the MAPQ 0 that depends on them is only a quarter
> of the total. The rest is ordinary multi-mapping on repeats, segmental
> duplications and the 2,385 decoys, which is physiological and not a defect.
>
> The right diagnostic metric is not global MAPQ 0 but MAPQ 0 **inside the MHC**,
> where the defect is concentrated. There the measured reduction is enormous, and
> the row has been added to the table.

If the ordering of the recalls straightens out, the diagnosis was right and is
demonstrated. If it does **not**, there was something else as well, and knowing
that is worth the phase anyway. Either way the result gets published.

**Cost.** About an hour, almost all of it `vcfeval`.

---

### Phase 5 — Re-deciding the hard filters

**What it covers.** A comparison of three scenarios on the new numbers: no
filters, current filters, filters without `SNP_MQ40`. The decision is taken on
measured F1, not by convention.

**Nothing needs re-running.** The ROC curves hap.py already produces
(`reports/giab/happy.roc.Locations.SNP.csv.gz` and the indel equivalent) contain
precision and recall at every possible threshold. The hard-filtered VCF keeps all
records and only annotates the `FILTER` column, so any subset can be derived after
the fact.

**Hypothesis to test — falsified on 30 July 2026.** I had written that, once the
upstream defect was removed, `SNP_MQ40` would stop biting on its own, because
MAPQ would no longer be artificially depressed. **That is wrong.** Measured on the
hard-filtered VCF of phase 3:

| Filter | Baseline | Iteration 2 | Delta |
|---|---|---|---|
| `SNP_MQ40` | 207,014 | 205,464 | **−0.7 %** |
| `SNP_SOR3` | 71,572 | 74,153 | +3.6 % |
| `SNP_QD2` | 46,312 | 50,428 | +8.9 % |
| `INDEL_QD2` | 4,448 | 4,954 | +11.4 % |
| PASS share | 94.38 % | 94.44 % | +0.06 % |

The faulty reasoning was assuming that the records marked by `SNP_MQ40` were the
same ones the ALT defect damaged. They were not.

Where the ALT defect struck — the MHC first of all — MAPQ collapsed to ~3, so
HaplotypeCaller discarded the reads *below its own threshold of 20* and **called
no variant at all**. Those sites never appeared in the VCF, and therefore could
not even be marked by a filter. They were absent, not filtered. After the fix
their mean MAPQ is ~58, so they get called **and** they clear `MQ ≥ 40` without
the filter touching them.

The 207,014 records marked by `SNP_MQ40` are a different population: sites in
ordinary repeats and segmental duplications, with intermediate MAPQ, where enough
reads cleared the threshold of 20 to allow a call but the root-mean-square MQ
stayed below 40. Those regions have nothing to do with ALT contigs, and indeed
the residual genome-wide MAPQ 0 is still 4.86 %.

**Consequence for this phase.** The question "do the SNP filters cost more than
they return?" stays **open and must be decided on hap.py**, not on these counts.
It is likely, though, that the answer has not changed, because the target
population of the dominant filter has not changed: in the baseline the filters
cost 41,350 true positives to remove 13,583 false ones.

A note of caution in the opposite direction: `SNP_QD2` (+8.9 %) and `INDEL_QD2`
(+11.4 %) grow more than total records (+2.91 %). QD is quality normalised by
depth, and more markings means more sites with marginal support — consistent with
the Ti/Tv of the added SNPs, around 1.95 instead of the 2.10 expected for
variants that are all true. Some of the gain will be false. How much, phase 4
says.

**Decided on 30 July 2026, on the phase 4 numbers.** The SNP hard filters still
cost more than they return, exactly as in the baseline:

| Type | F1 ALL | F1 PASS | Verdict |
|---|---|---|---|
| SNP | **0.9921** | 0.9883 | ALL wins: the filters remove 0.0038 of F1 |
| INDEL | 0.9924 | **0.9928** | PASS wins by a whisker: the filters are neutral |

Going from ALL to PASS on SNPs loses 40,409 true positives to remove 15,445 false
ones: **2.6 good ones for every error caught**, against 3.0 in the baseline.
Marginally less bad, but still a bad trade — and it is exactly what the analysis
of the mechanism predicted, given that the target population of the dominant
filter has not changed.

**Recommendation: turn off the SNP hard filters, keep the indel ones.**
Alternatively, replace them with thresholds derived from the ROC curves instead
of GATK's historical thresholds, which were tuned on data and depths different
from these.

**Cost.** No GPU compute.

---

### Phase 6 — Ploidy of chrX, chrY and chrM

**What it covers.** Re-calling HaplotypeCaller **only** over the haploid intervals
with `--ploidy 1`, excluding the pseudoautosomal regions which stay diploid, and
substituting those records in the final VCF. chrM has to be declared not
analysed: heteroplasmy requires a dedicated workflow, and presenting it diploid is
worse than omitting it.

**Cost.** About 30 minutes, because the intervals are small. It does not require
redoing the alignment.

**Acceptance criterion.** Zero heterozygous genotypes on chrY outside the PARs.

---

### Phase 7 — Documentation update

**What it covers.** `README.md` updated with two side-by-side columns, before and
after. The teaching value of the project is not in the final number: it is in the
complete cycle — a defect found with a benchmark, explained by a mechanism,
corrected, and the correction verified against a prediction stated in advance.

The limitations section needs updating too: the ones that were solved move into
the history, they do not get deleted.

---

## 4. Summary

| Phase | What | Cost | Blocks |
|---|---|---|---|
| 0 | Checkpoint provenance | code | everything else |
| 1 | No-ALT reference | ~4 GB of download | phase 2 |
| 2 | Smoke test | ~15 min | phase 3 |
| 3 | Full run | ~10 h 40 m | phase 4 |
| 4 | Benchmark and predictions | ~1 h | phases 5 and 7 |
| 5 | Hard filters, re-decided | none | phase 7 |
| 6 | X/Y/M ploidy | ~30 min | phase 7 |
| 7 | Documentation | writing | — |

Critical path: **phase 0 (half a day) → phases 1-2 (an hour and a half) → phase 3
(one night) → phases 4-5 (one morning)**. Phases 6 and 7 have no timing
constraints.

---

## 5. What we are not doing, and why

**DeepVariant.** Excluded by choice. It would not solve defect 1 anyway: the
damage happens in the aligner, and a variant caller only reads what the aligner
wrote. DeepVariant also filters on mapping quality (default
`min_mapping_quality = 5`), so it discards MAPQ 0 reads just the same. It would
help with defects 2 and 3, but it would introduce a new problem of benchmark
honesty: its WGS models are trained on HG001–HG007 excluding HG003, that is,
**HG002 is in the training set**. A genome-wide comparison on HG002 would be
rigged in its favour, and the only clean comparison would be restricted to
chr20–22, the only chromosomes always held out of training.

**Lowering HaplotypeCaller's MAPQ threshold** to recover the MHC without
realigning. It would call variants from reads that are genuinely in an uncertain
position: recall recovered and paid for in false positives. That is rigging the
number, not solving the problem.

**Adding the `.alt` file to the current reference.** In theory it is the correct
alt-aware route. In practice support in `pbrun fq2bam` changes between Parabricks
versions and `bwa-postalt.js` would have to run as a separate step outside
Parabricks. More fragile and less verifiable than the no-ALT reference, for the
same result.

**The pangenome** (`pangenome_germline` in Parabricks 4.7: Giraffe on GPU plus
pangenome-aware DeepVariant). It would attack the problem at the root, because
the reference contains the alternative variants instead of having to guess them.
But it requires the HPRC graph and more RAM than there is here. It stays a
hypothesis for an iteration 3, not for this one.

**Changing snpEff, known-sites and reference in the same run.** All three are
sensible improvements. Together they make the comparison unreadable.

---

## 6. Status

- [x] Phase 0 — checkpoint provenance · 29 July 2026
- [x] Phase 1 — no-ALT reference · 29 July 2026
- [x] Phase 2 — smoke test · 29 July 2026
- [x] Phase 3 — full run · 30 July 2026, 10 h 55 m, nine steps `PASS`
- [x] Phase 4 — benchmark and predictions · 30 July 2026, SNP F1 0.9921
- [x] Phase 5 — hard filters · 30 July 2026: off by default on SNPs in
      `scripts/postprocess.sh`, kept on indels, restorable with
      `HG002_SNP_HARD_FILTERS=on`
- [x] Phase 6 — X/Y ploidy · 30 July 2026: separate haploid call set for the
      regions outside the PARs, 110,099 genotypes, zero heterozygous. **chrM
      remains explicitly unaddressed**
- [x] Phase 7 — documentation · 30 July 2026: README with the two iterations
      side by side

---

## 7. Verified outcome of phases 0, 1 and 2

Run on 29 July 2026. Every number here comes from a file in the project.

### Phase 0

Prefix and volumes derive from the provenance key. Three tests passed:

| Test | Outcome |
|---|---|
| Reference with `_alt` and without `.alt` | the pipeline **refuses** to start |
| With `HG002_ALLOW_ALT_WITHOUT_ALT_FILE=1` | preflight passes, with a warning |
| Two identical invocations | same key, `97bd8ce8` |
| Change of reference | different key, `dba323b7` |

The plain-text manifest is written to `logs/runs/<id>/run_manifest.json`. The
baseline volume `hg002_work_v1` (91.1 GB) is intact.

Also added is the check that iteration 1 lacked: the preflight blocks a reference
carrying `_alt` contigs without the `.alt` file beside it.

### Phase 1

| Criterion | Expected | Measured |
|---|---|---|
| Total contigs | 2,580 | 2,580 |
| `_alt` contigs | 0 | 0 |
| HLA contigs | 0 | 0 |
| Decoys kept | 2,385 | 2,385 |
| Contigs added | 0 | 0 |
| NCBI MD5 | matching | FASTA and indexes OK |
| Regenerated `.fai` | identical to the published one | identical |
| `.dict` | 2,580 sequences | 2,580 |

786 contigs removed: 261 `_alt` plus the 525 HLA contigs that depend on them.
Names and lengths of chr1, chr6, chr10, chr20, chrX, chrY and chrM unchanged, so
dbSNP, snpEff and the GIAB BED stay valid without modification.

### Phase 2

Nine steps, thirteen executions counting the checks, all `PASS` in **5 m 54 s**
(`logs/runs/20260729_230546_smoke`). No checkpoint from the old reference was
reused: `2_fq2bam` really did run.

BAM header: 2,580 `@SQ`, zero `_alt`, zero HLA, 2,385 decoys. The HTML report
declares the reference actually used instead of a hand-written string.

**The main technical risk is behind us:** Parabricks 4.7 accepts the reference
with the `.fna` extension and the BWA indexes published by NCBI, without needing
to rebuild them.

### First evidence that the fix works

The smoke test uses the first million read pairs, spread across the whole genome,
so MAPQ is already comparable. Same input, the reference being the only
difference:

| Region | Reference | Reads | MQ0 | MQ0 % | Mean MAPQ |
|---|---|---|---|---|---|
| **MHC** `chr6:28.5–33.5 Mb` | with ALT | 538 | 492 | **91.45 %** | 3.4 |
| **MHC** `chr6:28.5–33.5 Mb` | no-ALT | **3,499** | 55 | **1.57 %** | **58.4** |
| whole chr20 (control) | with ALT | 45,475 | 2,558 | 5.63 % | 55.3 |
| whole chr20 (control) | no-ALT | 45,657 | 2,428 | 5.32 % | 55.5 |
| whole genome | with ALT | 1,994,995 | 129,569 | 6.49 % | — |
| whole genome | no-ALT | 1,994,892 | 96,931 | 4.86 % | — |

Three things to read in this table.

**Inside the MHC the mechanism is confirmed.** Mean MAPQ goes from 3.4 to 58.4
and MAPQ 0 reads from 91.45 % to 1.57 %. Above all, the primary locus receives
**6.5 times more reads** (538 → 3,499): with the old reference those reads were
smeared across the seven alternative MHC haplotypes instead of landing on the
real locus. It is exactly the causal chain described in chapter 1.

**On chr20 nothing changes**, as it should be. The fix acts where the defect was
and does not disturb the rest: the 12 variants called on `chr20:1-1,000,000` are
identical in the two runs, Ti/Tv included.

**Global MAPQ 0 barely falls** because ALT contigs cover a small fraction of the
genome. This falsified one of this plan's predictions, which was corrected in
phase 4 instead of being passed off as a hit.

Still unverified, and only verifiable after phase 3, are the metrics that require
full depth and comparison with the truth set: recall inside the MHC, genome-wide
recall, F1 and the ordering between SNPs and indels.

---

## 8. Phase 3 — what was observed during the run

Run `20260729_232130_full`, provenance key `53007e55`. Observations collected
while the pipeline was running, before the benchmark. **No number here is an
accuracy measurement**: that comes from phase 4.

### 8.1 Timings: the gap is where the reference changes the work

| Step | Baseline | Iteration 2 | Delta |
|---|---|---|---|
| `fq2bam` | 2 h 01 m | 2 h 13 m 05 s | **+9.8 %** |
| BQSR | 17 m | 18 m 29 s | **+8.7 %** |
| HaplotypeCaller | 4 h 41 m | 4 h 37 m 58 s | **−1.1 %** |
| **cumulative** | 6 h 59 m | 7 h 09 m | +2.5 % |

Three hypotheses were on the table for the initial gap: more useful work, the
GPU undervolt, memory pressure on the host. The last two predicted a **uniform**
slowdown, because they act on everything. The HaplotypeCaller figure clears them:
the gap is concentrated where removing the ALT contigs changes the amount of work
— alignment decides where every read goes, BQSR models more reads above threshold
— and disappears where it does not.

On HaplotypeCaller the two effects cancel out: 786 fewer contigs to traverse pull
downwards, while the previously blind regions that now generate real
ActiveRegions pull upwards.

### 8.2 Call set: +123,095 SNP alleles

| Metric | Baseline | Iteration 2 | Delta |
|---|---|---|---|
| Records | 4,953,729 | 5,098,062 | +144,333 (+2.91 %) |
| SNP alleles | 4,037,673 | 4,160,768 | **+123,095 (+3.05 %)** |
| Indel/complex alleles | 1,006,658 | 1,030,015 | +23,357 (+2.32 %) |
| Heterozygous | 3,104,375 | 3,212,050 | +107,675 (+3.47 %) |
| Homozygous alt | 1,849,354 | 1,886,012 | +36,658 (+1.98 %) |
| Ti/Tv | 1.928 | 1.929 | +0.001 |
| Mean DP | 40.57 | 40.88 | +0.31 |
| Records on `_alt` contigs | — | 0 | 0 expected |
| Contigs carrying variants | 1,913 | 1,689 | −224 |

The gain is of the same order as the deficit measured by the baseline benchmark
(101,963 true SNPs missing inside the high-confidence BED) and its composition is
the expected one: SNPs grow more than indels, because the ALT defect cost mostly
SNPs; heterozygous sites grow more than homozygous ones, because a heterozygous
site with one allele equal to the reference is invisible until there is coverage.

Inside the MHC window the SNPs called go from **301 to 28,766**. The number should
not be read as recall: it includes calls outside the high-confidence BED, and the
MHC is largely excluded from it precisely because it is hard.

**A cautionary signal.** Overall Ti/Tv does not move. Deriving the Ti/Tv of the
added SNPs alone from the totals gives about **1.95** — above the 1.93 of the
existing set, below the 2.10 of the truth set. It is a value derived from a
rounded ratio, so it carries an uncertainty of a few hundredths, but the direction
is clear: part of the gain is not real.

### 8.3 Hard filters: hypothesis falsified

See phase 5. In short: `SNP_MQ40` is immobile (−0.7 %), because the sites the ALT
defect damaged **were not filtered, they were absent**. `SNP_QD2` (+8.9 %) and
`INDEL_QD2` (+11.4 %) grow more than total records, a second confirmation that
some of the new calls have marginal support.

### 8.4 The most important result: the MHC has been recovered, but on a reference that cannot represent it

`PASS` records with predicted impact `HIGH` go from **613 to 762 (+149)**, that is
+24.3 % against +2.91 % of total records. It is not evenly distributed:

| Contig | Baseline | Iteration 2 | Delta |
|---|---|---|---|
| **chr6** | 16 | 64 | **+48** |
| chr19 | 40 | 65 | +25 |
| chr3 | 43 | 59 | +16 |
| chr11 | 42 | 57 | +15 |

**All 48 of the new chr6 records fall inside `chr6:28,510,120-33,480,577`, where
the baseline had none.** The MHC occupies 5 Mb, 0.16 % of the genome, and on its
own produces 32 % of the increase: an enrichment of roughly 200-fold.

The genes involved leave no doubt about what is happening:

| Gene | Records | | Gene | Records |
|---|---|---|---|---|
| HLA-DRB1 | 28 | | HLA-DQB1 | 6 |
| HLA-DRB5 | 14 | | HLA-B | 6 |
| MICA | 10 | | HLA-A | 6 |
| CCHCR1 | 9 | | TAP2 | 3 |
| TCF19 | 6 | | MUC21, OR12D1, PSORS1C1/2, HLA-J | 1–4 |

The predicted effects are 45 `frameshift_variant`, 12 `stop_gained`, 7
`splice_acceptor_variant`, 6 `splice_donor_variant`, 3 `stop_lost`. A pile of
loss-of-function calls in the classical HLA genes.

**And the technical quality is good:** mean DP 36.2, in line with the 35×
genome-wide, mean GQ 83.5, only 4 sites with DP < 10 and 3 with GQ < 30. The
variant caller is confident.

And that is precisely the point. **Technically solid, biologically meaningless as
losses of function.** HLA-DRB1, HLA-DRB5, HLA-A, HLA-B and HLA-DQB1 are the most
polymorphic genes in the human genome, with thousands of known alleles differing
from one another by dozens of substitutions and insertions. The GRCh38 primary
assembly carries **one arbitrary haplotype** of them. HG002's real HLA alleles
differ from that template so much that aligner and caller describe the difference
as a sequence of frameshifts and stop codons. The frameshift is not in HG002's
biology: it is in the reference, which is the wrong template.

The case of HLA-DRB5 is the most instructive: DRB5 exists only on some DR
haplotypes. Reads coming from DRB1, a paralogue identical for more than 90 %,
pile up on the reference locus of DRB5 and produce nonsense. Removing the ALT
contigs **favours** this effect, because those reads have nowhere else to go. The
`1/2` genotype on HLA-A at chr6:29,943,463 is a direct clue too: a multiallelic
site where both haplotypes of the sample differ from the reference.

**Three consequences.**

1. The fix did what it was supposed to: the MHC is no longer blind. That remains
   a success, and the phase 4 benchmark will quantify it.
2. The benchmark **will not penalise** these 48 records, because the MHC is
   largely outside the high-confidence BED. F1 can improve while those records
   stay artefacts. An improvement in F1 is not a promotion of the `PASS + HIGH`
   shortcut.
3. The `PASS + HIGH` list has to be read with one more caveat, now quantified:
   48 records out of 762, 6.3 %, sit in the MHC and are in all likelihood
   divergence from the reference, not loss of function. HLA needs dedicated
   typing tools, not variant calling on a linear reference, or a pangenome
   reference containing the alternative haplotypes instead of forcing them onto
   one.

It is the strongest argument in favour of the pangenome hypothesis for a possible
iteration 3, and it comes from a measurement, not from a preference.

### 8.5 Annotation: no improvement, as was not predicted

Records on contigs the snpEff database does not recognise go from 45,177 to
45,734. As a proportion of total records they barely fall, from 0.912 % to
0.897 %. I expected a clear drop, having removed 786 contigs: wrong, because
those records are dominated by the 2,385 decoys, which were **deliberately
kept**. The limitation stated in the README stands as it is.

### 8.6 Phase 4 — the benchmark

`hap.py` v0.3.12 with the `vcfeval` engine, NIST HG002 GRCh38 v4.2.1 truth set,
restricted to the high-confidence BED. Results in
`reports/giab/HG002_NovaSeq_40x_53007e55/`.

| Type | Filter | Metric | Baseline | Iteration 2 | Delta |
|---|---|---|---|---|---|
| SNP | ALL | Recall | 96.97 % | **99.26 %** | **+2.29 pp** |
| SNP | ALL | Precision | 99.33 % | 99.16 % | −0.18 pp |
| SNP | ALL | **F1** | 0.9814 | **0.9921** | **+0.0107** |
| SNP | ALL | False negatives | 101,963 | **24,740** | **−75.7 %** |
| SNP | ALL | False positives | 21,931 | 28,440 | +29.7 % |
| INDEL | ALL | Recall | 97.36 % | **99.19 %** | **+1.84 pp** |
| INDEL | ALL | Precision | 99.40 % | 99.29 % | −0.11 pp |
| INDEL | ALL | **F1** | 0.9837 | **0.9924** | **+0.0088** |
| INDEL | ALL | False negatives | 13,895 | **4,238** | **−69.5 %** |

**The ordering of the recalls has straightened out.** In the baseline indels
(97.36 %) beat SNPs (96.97 %), which is the opposite of how it should be. Now SNPs
are at 99.26 % against 99.19 % for indels: SNPs are the easy case again. This was
the plan's main prediction and it is confirmed.

**The MHC.** SNP recall inside `chr6:28,510,120-33,480,577`:

| | True | Recovered | Lost | Recall |
|---|---|---|---|---|
| Baseline | 20,177 | 301 | 19,876 | **1.49 %** |
| Iteration 2 | 20,177 | 19,658 | 519 | **97.43 %** |

**False negatives per chromosome**, the four hardest hit by the defect:

| Contig | Baseline | Iteration 2 |
|---|---|---|
| chr6 | 23,093 (10.39 %) | 1,425 (**0.64 %**) |
| chr17 | 7,585 (8.80 %) | 630 (**0.73 %**) |
| chr15 | 6,657 (6.47 %) | 1,543 (1.50 %) |
| chr19 | 3,714 (5.26 %) | 394 (**0.56 %**) |
| **genome** | 101,963 (3.03 %) | **24,740 (0.74 %)** |

chr6 goes from worst chromosome in the genome to one of the best. The improvement
is not confined to the ALT-rich chromosomes: chr1 also falls from 1.64 % to
0.91 % and chr2 from 1.34 % to 0.70 %, because alternative contigs exist a little
everywhere.

### 8.7 The predictions, verified one by one

| Prediction | Expected | Measured | Outcome |
|---|---|---|---|
| MAPQ 0 genome-wide | ~4.9 % (corrected) | 4.86 % | **right** |
| MAPQ 0 inside the MHC | 1–2 % | 1.57 % | **right** |
| SNP recall inside the MHC | > 95 % | 97.43 % | **right** |
| SNPs lost on chr6 | ~2 % | 0.64 % | **right**, better in fact |
| SNP recall genome-wide | ~99 % | 99.26 % | **right** |
| Ordering SNP > indel | yes | yes | **right** |
| Part of the gain is false | yes | +6,509 FP | **right** |
| SNP F1 (ALL) | 0.993–0.995 | 0.9921 | **wrong**, optimistic |
| `SNP_MQ40` stops biting | yes | −0.7 % | **wrong** (see phase 5) |

Seven out of nine. The two wrong ones stay written down, with the reason.

**The cost in precision, quantified.** The gain brings 77,223 more true positives
and 6,509 more false positives: **11.9 true for every false one**. The two
cautionary signals collected during the run — Ti/Tv of the added SNPs around 1.95
and `QD2` filters on the rise — correctly indicated that part of the gain was
false, and the benchmark says how much: not a lot.

**What remains.** F1 0.9921 against the ~0.995 a GATK pipeline should reach at
this depth. The bulk of the gap is closed, not all of it. Still on the table are
ordinary segmental duplications, the residual 4.86 % of MAPQ 0, and a BQSR
known-sites set limited to dbSNP alone.
