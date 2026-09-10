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
