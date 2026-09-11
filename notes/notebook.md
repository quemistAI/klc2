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

## Part 1.4 CORRECTED — hg38 coordinates via BLAT

First liftover attempt gave chr11:66,257,086-66,257,302. **Wrong.**
Sequence at those coords did not match the confirmed hg19 217-mer.
Its revcomp mapped to hg19 66,024,426-66,024,642 — 131 bp off.

**BLAT of the confirmed hg19 217-mer against hg38:**
chr11:66,256,955-66,257,171, MINUS strand, 100% identity,
full-length (query 1-217), span 217. Unique hit; next best is 39 bp.

**-> Correct hg38 interval: chr11:66,256,955-66,257,171**
BED: chr11  66256954  66257171

**The locus is INVERTED between GRCh37 and GRCh38.** hg19 plus
strand = hg38 minus strand. Verified: under this inversion,
hg19 66,024,426-66,024,642 maps to hg38 66,257,086-66,257,302,
which is exactly what the bad liftover returned.

**Rule for this project:** canonical sequence = hg19 plus-strand
217-mer (raw/klc2_del_217bp.fa). Matches published amplicon
orientation. Do not substitute the hg38 plus-strand version.
Verify any coordinate conversion by BLAT, not liftover alone.

**hg38 sequence verified:** chr11:66,256,955-66,257,171 plus strand
is the exact revcomp of the validated hg19 217-mer (217/217, no
mismatches). Confirms BLAT coords and the inter-build inversion.
Saved as raw/klc2_del_217bp_hg38.fa — this is the FIMO input.
GC content 54.4% (genome avg ~41%). Check CpG island overlap in Part 3.

**KLC2-AS1**: chr11:66,264,777-66,265,666 (hg38), ENST00000530805.1
Distance from my interval (ends 66,257,171): ~7.6 kb.
Too distant for a promoter-proximal antisense model. Demoted to
footnote — would need Hi-C/4C contact evidence to revive.
**KLC2-AS2**: [to check]

### Question 3 — CpG island: NO OVERLAP (but shore-adjacent)

**CpG island (hg38):** chr11:66,257,440-66,258,782
- Size 1,343 bp | CpG count 120 | 66.1% C+G | obs/exp 0.83
- Strong island; likely the KLC2 promoter (confirm strand/TSS)

**My interval:** chr11:66,256,955-66,257,171
**Gap:** 268 bp (66,257,172-66,257,439). No overlap.

**Interpretation — element sits in the CpG island SHORE.**
Shores (~2 kb flanking an island) are where most tissue-specific
differential methylation occurs; island cores are typically
unmethylated across tissues. A shore location is consistent with
the lineage-restricted activity implied by Melo 2015 (blood vs
fibroblast/MN difference).

GC contrast: element 54.4% vs island 66.1% — element is GC-rich
but distinct from the island core. Consistent with shore.

**Follow-up:**
- Confirm KLC2 strand and TSS position; is this island the promoter?
- Distance from element to KLC2 TSS?
- Check Roadmap/ENCODE WGBS methylation over the shore in neural
  vs blood in Part 4 — shore methylation difference would be a
  direct test.
- Revises earlier note: element is NOT in a CpG island. Do not
  weight CxxC-domain factors heavily in motif candidate ranking.

### Part 3 — MAJOR reframe: bidirectional promoter architecture

**ENSG00000255320 (ENST00000791493.1), GENCODE V50:**
chr11:66,244,570-66,257,625 hg38 | minus strand | 3 exons | 13,056 bp
| **No protein** -> lncRNA. TSS = 66,257,625 (high coord, minus strand).
Likely = KLC2-AS2. **My element lies INSIDE this transcript (intronic).**

**Locus architecture (L->R):**
- element      66,256,955-66,257,171
- KLC2 TSS     ~66,257,270 (PLUS strand, transcribes right)
- CpG island   66,257,440-66,258,782
- lncRNA TSS   66,257,625 (MINUS strand, transcribes left)

-> **Divergent/bidirectional promoter.** Element is ~100 bp upstream
of KLC2 TSS and ~450 bp from the lncRNA TSS.

**CTCF (ReMap filtered):** peaks present in window but concentrated
over the KLC2 promoter / CpG island, right of the element. Only wide
peak tails reach into the interval. Element is NOT a dedicated CTCF
site. Insulation-loss model demoted.

**NEW LEADING HYPOTHESIS (revises Hypothesis E):**
Element is an enhancer for the antisense lncRNA. Antisense
transcription represses KLC2 in cis. Deletion -> less lncRNA ->
more KLC2. Reconciles the S1B "enhancer" ChromHMM call with the
observed overexpression, which no prior work has explained.
DIRECT TEST: measure ENSG00000255320 in patient fibroblasts.

**Caveat:** GENCODE track says "14 items filtered out" — confirm no
transcripts are hidden before finalizing this architecture.

**GENCODE filter cleared (Tagged Sets: All, pseudogenes on):**
13 lncRNA isoforms now visible, ALL ENSG00000255320, all minus
strand. TSS cluster spans ~66,257,400-66,257,625 (multiple
alternative starts, not one). **No transcript initiates inside
the element.** Bidirectional-promoter model holds.

13 annotated isoforms = well-supported locus, not a spurious call.
Element is intronic to all of them.

KLC2 side: only 2 isoforms have extended 5' ends reaching left
toward the element; most start ~66,257,600+. State element-to-TSS
distance per isoform, not as a single figure.

Conservation unchanged: flat phyloP, mouse/dog/elephant alignment
absent across the interval. Primate-specific.

**Antisense identity resolved — THREE distinct loci:**
| ENSG00000255320 | 66,244,570-66,257,625 | CONTAINS element | no symbol |
| KLC2-AS1        | 66,264,777-66,265,666 | ~7.6 kb away     |
| KLC2-AS2        | 66,267,635-66,268,129 | ~10.5 kb away    |
(KLC2-AS2 = ENST00000533287.1, uc058dut.1)

**ENSG00000255320 is NOT KLC2-AS2.** No HGNC symbol assigned.
**Use the ENSG ID consistently. Never write KLC2-AS2 for it.**
Unnamed 13-kb, 13-isoform lncRNA over a disease gene promoter —
suggests an under-studied locus.

### Step 3 — GTEx ENSG00000255320: HYPOTHESIS E REFUTED

Expressed broadly ~2-5 TPM across ~50 tissues. Exceptions:
- **ALL brain regions LOWEST on the plot (~0-1 TPM)** — cerebellum,
  cortex, hippocampus, hypothalamus, putamen, substantia nigra,
  amygdala, spinal cord
- Cultured fibroblasts ~4 TPM (high end)
- Whole blood ~2-3 TPM (ordinary)

**Refutes antisense cis-repression model on two counts:**
1. lncRNA is LOWEST in brain, but Melo found KLC2 OVEREXPRESSED in
   motor neurons. Wrong direction.
2. Fibroblasts ≈ blood here, but Melo's central finding is that
   those tissues DIFFER for KLC2.

**Hypothesis E dropped. Do not revive.**

Reverting to primary model: element = lineage-restricted silencer
acting on KLC2 directly. Unaffected by today's findings.
Bidirectional-promoter architecture stands as a location finding,
not a mechanism.

Unexplained side observation (note only, do not build on):
brain-specific depletion of this lncRNA at a locus where a
non-coding deletion causes a brain disease.

### PREDICTION for Part 4 — written before running the analysis

Model: the 216-bp element is a silencer with lineage-restricted
activity, acting on KLC2 directly.

I predict the Roadmap ChromHMM overlap will show:
- Neural (E073, E081, E053, E054): [E073]
- Fibroblast (E055, E056): [either, possible E056]
- Blood (E062, E116): [E116]
Melo's data imply fibroblasts should pattern with neural, not blood.
If fibroblasts pattern with blood instead, the 2015 fibroblast
overexpression result needs another explanation.

Confidence: [medium, fibroblasts and blood may be closer than with neural according to the GTEx tissue analysis]

Caveat I already expect: motor neurons are absent from Roadmap.
A null or quiescent result may reflect missing cell types rather
than absence of function.
