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
nano notes/notebook.mdPublished deletion (66,024,557-66,024,773) sits inside this amplicon
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

### Roadmap EID metadata verification — 2026-09-11

```
E053	#FFD924	Cortex derived primary cultured neurospheres
E054	#FFD924	Ganglion Eminence derived primary cultured neurospheres
E055	#FF9D0C	Foreskin Fibroblast Primary Cells skin01
E056	#FF9D0C	Foreskin Fibroblast Primary Cells skin02
E062	#55A354	Primary mononuclear cells from peripheral blood
E073	#C5912B	Brain_Dorsolateral_Prefrontal_Cortex
E081	#C5912B	Fetal Brain Male
E116	#000000	GM12878 Lymphoblastoid Cells
```
## Part 4 RESULT

| EID | Tissue | Group | State |
|-----|--------|-------|-------|
| E053 | Neurosphere ctx | Neural | 2_TssAFlnk |
| E054 | Neurosphere GE | Neural | 1_TssA |
| E073 | DLPFC | Neural | 2_TssAFlnk |
| E081 | Fetal brain | Neural | 2_TssAFlnk |
| E055 | Fibroblast | Fibro | 2_TssAFlnk |
| E056 | Fibroblast | Fibro | 2_TssAFlnk + 1_TssA |
| E062 | PBMC primary | Blood | **11_BivFlnk** |
| E116 | GM12878 | Blood | 2_TssAFlnk + 1_TssA |

**Prediction WRONG.** Predicted repressive states in neural/fibro.
Got active promoter states in 7/8. No ReprPC anywhere.

**Single lineage difference: E062 primary blood = bivalent (poised).**
Aligns with Melo: active promoter in tissues where deletion raises
KLC2 (neural, fibroblast); poised where it doesn't (blood).
Different mechanism than hypothesized, same lineage split.

**Caveat 1:** E116 (also blood) shows active states, contradicting
E062. GM12878 is EBV-transformed — chromatin altered by
transformation. E062 (primary) more trustworthy, but n=1.

**Caveat 2 — RESOLUTION.** Segments are 200-1800 bp; element is
217 bp sitting ~100 bp from the KLC2 TSS. ChromHMM 200-bp bins
cannot separate element from promoter here. These calls largely
report PROMOTER state. Must state this explicitly in any writeup.

**Next:** (1) expand to 6-8 primary blood epigenomes (E029-E051) to
test whether bivalency reproduces. (2) Move to continuous H3K27me3 /
H3K4me1 / H3K4me3 bigWig signal for base-level resolution.

## Part 4 FINAL — expanded panel (n=16)

**Threshold (prespecified): 6+/8 blood bivalent. RESULT: 4/9. FAILED.**

Myeloid:  E029 EnhBiv | E030 ReprPC+EnhBiv | E124 EnhBiv+BivFlnk  -> 3/3
Mixed:    E062 PBMC BivFlnk (contains monocytes)
Lymphoid: E032, E034, E046, E047 all TssAFlnk/TssA              -> 0/4
Stem:     E035 TssAFlnk
Neural:   E053 E054 E073 E081 all active
Fibro:    E055 E056 all active
GM12878:  E116 active (transformed)

**Unhypothesized pattern: myeloid = bivalent/Polycomb, lymphoid =
active.** E124 independently replicates E029. HSC active, so not a
stem effect. Coherent, but POST HOC — not a finding until
prespecified and tested in an independent panel.

**Does NOT explain SPOAN.** Melo: neural+fibro overexpress, blood
doesn't. Here neural, fibro, and lymphoid blood are all
indistinguishable (active promoter). The lineage boundary found runs
WITHIN blood — orthogonal to the disease.

**Technical:** E030 ReprPC segment (66,022,200-66,024,600) overlaps
only ~44 bp of the element before switching to EnhBiv. Score as
bivalent.

**Conclusion: ChromHMM at 200-bp resolution cannot answer this
question** at an element ~100 bp from an active TSS. Segment calls
are

### Part 5 threshold (prespecified)

The element is distinguishable from the promoter if H3K4me3 signal
over the 217 bp is <50% of its peak value at the KLC2 TSS in the
same epigenome.

The element has its own regulatory signature if H3K4me1 exceeds
H3K4me3 over the interval in at least 3 of 5 epigenomes.

If H3K4me3 is uniformly high across element and promoter in all
five, the element is not separable from the promoter by any
reference chromatin data, and the reporter assay becomes the only
route.

### Part 5b prediction — written before downloading

Testing whether the myeloid H3K27me3 signal replicates.

Predicted mean fold-change over the 217-bp element:
- Myeloid (E030 neutrophil, E124 CD14+ monocyte):
    H3K27me3 ____  H3K4me3 ____
- Lymphoid (E032 B cell, E034 T cell, E047 CD8+ naive):
    H3K27me3 ____  H3K4me3 ____

Threshold for calling it replicated:
H3K27me3 > 10 in all 3 myeloid AND < 3 in all lymphoid.
H3K4me3 < 3 in all myeloid AND > 8 in all lymphoid.

If lymphoid looks like myeloid, the boundary is blood-vs-other,
not myeloid-vs-lymphoid, and E029 was not special.

## Part 5b RESULT — myeloid H3K27me3 REPLICATES

H3K27me3 mean fold-change over the 217-bp element:
  Myeloid:    E029 19.45 | E030 29.44 | E124 35.05
  Lymphoid:   E032 0.14  | E034 0.81  | E047 0.44
  Fibroblast: E055 1.12
  Neural:     E073 1.28

**No overlap. ~15x gap between lowest myeloid and highest other.**
Threshold (>10 all myeloid, <3 all lymphoid): MET.

H3K4me3 criterion (<3 in myeloid): **FAILED** — E030 7.26,
E124 10.58. High H3K27me3 + substantial H3K4me3 = BIVALENT, not
silenced. Matches ChromHMM 12_EnhBiv / 11_BivFlnk independently.

**Refined claim: the element carries a bivalent chromatin domain
specific to the myeloid lineage.** Confirmed in two independent
data types (segmented states, continuous signal).

Caveat: E047 low across all marks — possible weak-signal epigenome.
Caveat: still POST HOC and still orthogonal to SPOAN — neural and
fibroblast, the tissues where deletion raises KLC2, are
indistinguishable from lymphoid here.

## Part 5b SPATIAL — broad domain, not a peak

H3K27me3 in myeloid (E029, E030, E124) forms a BROAD ~1 kb domain
from ~66,023,800 rising through the element and terminating
sharply at ~66,024,800 where H3K4me3 promoter signal begins.
**Not a discrete peak over the 217 bp.**

**Claim must be narrowed to:** the element lies within a
myeloid-specific H3K27me3 domain upstream of the KLC2 promoter.
Any sequence in that ~1 kb would show the same signal; the 217 bp
is NOT distinguished within it.

H3K4me1 in myeloid (6-7) ≈ fibroblast (6.0). The "bivalent" reading
rests on H3K4me3, which is visibly promoter signal extending
leftward. Soften bivalency language.

Supports Part 3 architecture: sharp H3K4me3 onset at the right edge
of the element in all 8 epigenomes = the promoter boundary.

E047 weak epigenome (y-max 30 vs 125 elsewhere). Downweight.

### CTCF motif — ReMap binding evidence: SUPPORTED

Motif (FIMO pos 81-111) = hg38 chr11:66,257,035-66,257,065,
center 66,257,050.

| Dataset | Position | Width | Midpoint | Dist |
|---|---|---|---|---|
| GSE117508 cortical-interneuron KCl-neg | 66,256,969-66,257,149 | 181 | 66,257,059 | 9 |
| GSE117508 cortical-interneuron KCl-pos | 66,256,964-66,257,090 | 127 | 66,257,027 | 23 |
| ENCSR000AQU DND41 | 66,256,926-66,257,247 | 322 | 66,257,087 | 37 |
| GSE116862 hESC_D5 | 66,256,971-66,257,489 | 519 | 66,257,230 | 180 |
| ENCSR072EUE OCI-Ly1 | 66,257,014-66,257,503 | 490 | 66,257,259 | 209 |

**The two NARROWEST peaks are centered on the motif, and both are
NEURAL (cortical interneuron).** Broad peaks (hESC, OCI-Ly1) are
centered on the KLC2 promoter instead. Narrow peak = precise
localization.

-> FIMO's top hit has independent ChIP-seq support in neural cells.
-> **Insulation-loss model reinstated as best-supported mechanism.**

CAVEATS: both cortical-interneuron datasets are from one study
(GSE117508, KCl +/-), so one experiment observed twice, not two
independent replications. CTCF binding =/= insulation; demonstrating
that requires Hi-C/4C contact data, out of scope here.

NOTE: this came from neural ChIP-seq — exactly the data the 2015
analysis lacked.

| TF | FIMO p | ReMap status |
|----|--------|--------------|
| CTCF | 1.6e-5 | **SUPPORTED** — 2 neural peaks centered 9 and 23 bp from motif |
| REST | 8.4e-5 | **NOT supported** — 8 datasets in window (incl. neural, hippocampus); peaks cross the element with their edges, none centered on the motif at 66,256,979-998 |
| ZNF701 | 1.6e-5 | Not assessed — minimal ChIP-seq coverage in ReMap |
| ZKSCAN3 | 4.4e-5 | Not assessed — same |
| ZNF528 | 1.8e-5 | Not assessed — same |

CTCF is the only FIMO candidate with binding evidence localized to
its motif. The REST negative is informative rather than absent:
REST IS bound in this window in neural tissue, just not at the
predicted site.

## 2026-09-14 — FRAMING CORRECTION after review by U. Melo (first author, Melo et al. 2015)

**Source:** email correspondence, 2026-09-14.

### The error

I had framed the gap as: "Supplementary Fig S1B calls the element a
probable enhancer, but deleting it RAISES KLC2 expression, so the
enhancer call contradicts the direction of effect."

**This is an error.** ChromHMM states and ENCODE cCRE classes
describe CHROMATIN STATE, not direction of transcriptional effect.
A repressive sequence sitting inside an active promoter reads as
active chromatin because the promoter is active. No contradiction.

**The direction was never in question.** Melo et al. already showed
the element is repressive: reporter constructs LACKING the 216 bp
gave higher activity in HEK293T, U87MG and motor neurons, and
patient fibroblasts and iPSC-MNs overexpress KLC2.

### Corrected framing

The 216-bp element is an ESTABLISHED repressive sequence within the
KLC2 promoter. What is unknown is WHICH SEQUENCE FEATURE within the
216 bp mediates that repression. No dataset can
resolve a 216-bp element ~100 bp upstream of the TSS of a broadly
expressed gene.

### New research question

Which sequence feature within the 216-bp SPOAN element mediates its
repressive effect on KLC2, and can that feature be identified from
base-resolution data?

### Consequences

- Fig 3 result is EXPECTED, not a failed prediction. An element in an
  active promoter of a broadly expressed gene should read active in
  every tissue. Make into a resolution result.
- Motif scanning should target REPRESSORS (REST/RE1, KRAB-ZFPs), not
  enhancer-associated factors. My original priority list was right, but for the wrong reason.
- TSS distance: ~100 bp upstream of the nearest annotated TSS,
  ~600 bp upstream of the MANE TSS.

### Reviewer's other points (see revision plan)

1. CTCF p=1.6e-5 for a 19-bp motif is modest. Peak width reflects
   experiment, not specificity. Must check whether the site is
   constitutive (GM12878, K562, HepG2) and whether cohesin
   (RAD21, SMC3) also occupies.
2. Cortical interneurons are NOT the target cell. SPOAN is motor
   neuron / peripheral nerve / optic. Use tibial nerve, spinal cord,
   GTEx data- Fibroblasts are also a valid primary model.
3. Check the DELETION JUNCTION — a deletion can create a motif as
   well as remove one.
4. Reporter design: WT / 216-bp deletion / CTCF core mutated (4-5 bp).

## 2026-09-14 — Part 2 & 3: cohesin negative, junction negative

### Cohesin co-occupancy at the CTCF motif: NEGATIVE

CTCF motif = chr11:66,257,035-66,257,065 (hg38), centre 66,257,050.

**SMC3** (ReMap, hg38, window chr11:66,256,455-66,257,671): 2 peaks.
  - neural: chr11:66,257,081-66,257,544. Starts **16 bp past the 3' end
    of the CTCF motif**; midpoint 66,257,313, i.e. **262 bp from the
    motif centre**, over the KLC2 promoter.
  - peripheral-blood-neutrophil: starts 66,257,395. No overlap.

**RAD21**: 8 peaks (neural, lymphoblast x2, HEK293 x2, HeLa-S3, HAP1,
SLK). Leftmost (neural) starts 66,257,158, at the element boundary.
None centred inside the element.

**Conclusion:** no cohesin peak overlaps the CTCF motif. Cohesin is
positioned over the promoter. CTCF sites doing architectural work are cohesin co-occupied, so this argues the element's CTCF
site is **not** a loop anchor.

### Consequence: insulation-loss model demoted (third revision)

Track record of this hypothesis:
  1. Part 3 — demoted: ReMap CTCF peaks looked promoter-centred
  2. FIMO — revived: CTCF was the top-ranked motif (p=1.6e-5)
  3. Cohesin — demoted again: no RAD21/SMC3 at the motif

**Aim 2 (mutate the CTCF core) remains justified**, but the
justification changes. NOT "the deletion removes a boundary."
Instead: CTCF is the highest-ranked motif in the element, making it the best single candidate
for the repressive feature.

### CTCF occupancy breadth

Full CTCF track, ~78 peaks in the window. **Only 8 reach the element**:
heart, DOHH2, DND41, OCI-Ly3, cortical-interneuron x2, hESC, OCI-Ly1.
The other ~70 start right of the element and cluster on the promoter
and CpG island.

GM12878 and HepG2: **absent from the window entirely**. K562: two
peaks, both right of the element.
-> **Site is NOT constitutive.**

BUT: of the 8 element-overlapping peaks, only 2 are neural. Three are
B-cell lymphoma lines (DOHH2, OCI-Ly3, OCI-Ly1). And neural cell types
are well represented in the promoter-only group (SH-SY5Y, neural,
nerve, retina, anterior-temporal-cortex).
-> **Correct claim: cell-type-restricted, NOT neural-specific.**
-> Peak width is not evidence of specificity.

### Deletion junction motif scan: NEGATIVE

Built 200-bp junction (100 bp each flank, element removed) ->
raw/klc2_junction_200bp.fa. Assertion confirmed no element sequence
leaked in. FIMO vs JASPAR2024 CORE vertebrates, p<1e-4.

**No motif spans the junction point.** The deletion does not create a
new transcription factor binding site at the breakpoint.

Rules out the gain-of-motif mechanism. The
deletion appears to act by REMOVING a feature, not creating one.

### Still to do
- Re-scan the 216 bp for repressors specifically (REST/RE1, KRAB-ZFPs)
- q-value for the CTCF hit
- DNase/DHS footprints over the element
- GTEx tibial nerve + spinal cord; search for iPSC-MN ATAC/ChIP
- Revise Fig 5 (drop neural colour coding and width-based sort)

### Hi-C (3D Genome Browser, hg38, 5 kb, GM12878 + IMR-90): NEGATIVE

Window chr11:66,150,000-66,400,000. Three loops called:
| anchor 1 | anchor 2 | score | cell |
| 66,060,001-66,070,000 | 66,230,001-66,240,000 | 229 | GM12878 |
| 66,060,001-66,070,000 | 66,280,001-66,290,000 | 192 | GM12878 |
| 66,080,001-66,090,000 | 66,180,001-66,190,000 |  26 | IMR-90  |

**Element (66,256,955-66,257,171) is in none of them.** Sits in the
~40 kb gap between the 66,230-240 kb and 66,280-290 kb anchors:
17 kb past one, 23 kb before the next.

Boundary position differs between cell types; tens of kb from
the element in both.

IMR-90 = lung fibroblast, the most disease-relevant cell type
available (patient fibroblasts overexpress KLC2). Only one weak
loop (score 26), neither anchor near the element.

**LIMITATION: 5 kb resolution, 10 kb loop anchors, vs a 217-bp
element.** Hi-C cannot address this element directly. It can say the
LOCUS is not a loop anchor. Consistent with the cohesin negative.

-> Insulation-loss model now has three independent negatives:
   no cohesin at the motif, no loop anchor at the locus, and CTCF
   bound in only 8/78 cell types at the element.
 
### DNase + H3K27ac (UCSC layered, averaged by organ/tissue, hg38)

**DNase:** near baseline across the left two-thirds of the element.
Rises from ~66,257,120, peaks ~66,257,300 over the promoter.
-> Element is NOT in accessible chromatin; the adjacent promoter is.
-> Further argues against broad CTCF occupancy at 66,257,035-66,257,065
   (CTCF binding generally requires accessibility). Consistent with
   8/78 cell types.

**H3K27ac:** rises approaching the element, DIPS within it
(~66,257,050-66,257,200), rises again to peak past the right edge.
-> **First base-resolution feature localized to the element rather
   than the promoter.** A local trough in active marking inside an
   otherwise active promoter — the expected signature of a repressive
   sequence embedded in a promoter.
-> Weak: averaged data, modest amplitude, not quantified. Treat as an
   observation, not a result.

Caveat: both tracks are organ/tissue averages, not single cell types.
No footprint-level (base-resolution protection) data located.

### GTEx KLC2 (median TPM) — REINTERPRETS the Melo blood result

Nerve - Tibial            24.43  (n=670)
Brain - Spinal cord c-1   19.53  (n=204)
Cells - Cultured fibro    19.03  (n=652)
**Whole Blood              2.77  (n=803)**
Cerebellum / cereb hemi   ~190   (highest in the body)

**Blood is 7-9x lower than every disease-relevant tissue.**
-> Melo's "no change in blood" is likely an EXPRESSION-LEVEL effect,
   not a chromatin or element-activity difference: KLC2 is barely
   expressed in blood, so there is little baseline for a repressive
   element to act.
-> Explains why Fig 3 found no chromatin difference along the disease
   axis. The neural-vs-blood contrast was never a chromatin question.

Note: cerebellum ~190 TPM, ~10x the affected tissues. KLC2 is highest
where SPOAN spares.

### H3K27ac over element vs 200-bp flanks (Roadmap fold-change, hg19)
E073 neural     element 5.54 | flank 4.65 | ratio 1.19  (no dip)
E055 fibroblast element 5.52 | flank 7.99 | ratio 0.69  (dip)
E029 monocyte   element 1.81 | flank 3.01 | ratio 0.60  (dip)

Dip in fibroblast and monocyte, absent in prefrontal cortex.
CAVEAT: flanks are asymmetric (E055 left 1.49 vs right 14.48) because
the promoter abuts the right side. Ratio partly reflects promoter
proximity, not a symmetric trough. Weak observation.

### Motor neuron accessibility
ENCODE: 10 human motor neuron ATAC-seq experiments (Snyder lab, all
released, ALS donor-derived). 0 DNase-seq. Coverage of the element
not yet checked.

