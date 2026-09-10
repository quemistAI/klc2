# KLC2 / SPOAN Project Notebook

## 2026-09-10 — Environment setup

**Goal:** get the computational environment running.

**Done:**
- Installed Miniconda, created env `klc2` (python 3.11)
- Added channels: conda-forge > bioconda > defaults, strict priority
- Installed: bedtools, pandas, numpy, matplotlib, seaborn,
  biopython, pybedtools, pyBigWig, pyGenomeTracks, jupyterlab
- Folder structure created, git repo initialized and pushed
  to github.com/quemistAI/klc2 (private)

**Versions:** see notes/environment.txt

**Confused about:** nothing yet.

**Next:** Part 1 — resolve the 216 vs 217 bp discrepancy against
Melo 2015 supplementary, then liftover hg19 -> hg38.

## 2026-09-10 — Part 1: starting interval

**Goal:** get the published SPOAN deletion into a BED file.

**Source:** Melo et al. 2015, HMG 24(24):6877-6885.
Published interval, GRCh37: chr11:66,024,557-66,024,773

**Ran:**
printf "chr11\t66024556\t66024773\tKLC2_SPOAN_del\n" > raw/klc2_del_hg19.bed
cat -A raw/klc2_del_hg19.bed

**Got:** chr11^I66024556^I66024773^IKLC2_SPOAN_del$
Tabs confirmed.

**Note:** start is 66024556 not ...557 — BED is 0-based, so start
minus one, end unchanged.

**Open question:** end - start = 217, but the deletion is called
216 bp everywhere, and ClinVar's c.-451_-235del is also 217.
Off-by-one somewhere. Need to check Melo supplementary for the
actual breakpoint sequence before designing any constructs.

**Next:** Part 1.3 — resolve the 216/217 question, then liftover.

### In-Silico PCR result (UCSC, hg19)

Primers from Table S3 pair 1 → **chr11:66,024,368-66,024,831, 464 bp**
Matches the wild-type band in Fig S3C exactly. Primers confirmed.
Published deletion (66,024,557-66,024,773) sits inside this amplicon
with 189 bp of 5' flank and 58 bp of 3' flank.

**Breakpoint analysis (mapped published coords onto amplicon seq):**
- 5' context 66,024,547-66,024,566: ACTACGTTTGGCTCTGAGAC
- 3' context 66,024,764-66,024,783: GCAGTTTCCTGACCCTTCCC
- First base of deletion (66,024,557) = G
- Base immediately after deletion (66,024,774) = G
  → **1-bp microhomology.** 66,024,557-66,024,773 and
    66,024,558-66,024,774 delete identical sequence. Position is
    ambiguous by one base (standard HGVS 3'-shift situation).

**But microhomology does not explain the size discrepancy** —
shifting changes position, not length. Both versions = 217 bp.
Gel-derived size = 216 bp. So there is a genuine off-by-one in the
published coordinates; one endpoint is likely 66,024,772 (or start
66,024,558).

**Cannot resolve definitively without patient junction sequence.**
Not in the supplement. Options: contact corresponding author
(Melo / Santos), or Sanger a patient sample if ever available.

**Decision:** proceed with the published 217-bp interval as the
analysis window. Immaterial for chromatin state, conservation, and
motif scanning. MUST be resolved before reporter construct design.

## Part 1.4 — LiftOver hg19 -> hg38

chr11:66,024,557-66,024,773 (GRCh37)  -> raw/klc2_del_hg19.bed
chr11:66,257,086-66,257,302 (GRCh38)  -> raw/klc2_del_hg38.bed
Span preserved: 217 bp both builds. Single contiguous block,
no split. Offset ~233 kb.

Reminder: Roadmap ChromHMM work stays in hg19. Everything else
(ENCODE, UniBind, ReMap, phyloP, JASPAR) uses hg38.
