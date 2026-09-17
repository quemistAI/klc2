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

## 2026-09-14 — Motor neuron ATAC-seq at the element

**Source:** ENCODE, human motor neuron ATAC-seq, Snyder lab (Stanford),
all released. 10 experiments available; 0 DNase-seq in motor neuron.
**All donors have amyotrophic lateral sclerosis.** No healthy-control
motor neuron ATAC identified.

Files viewed: ENCFF345PTN (rep 1,2), ENCFF638TCZ (rep 1,2),
ENCFF576QGF (rep 1), ENCFF721NVE (rep 1), ENCFF794PRL (rep 2),
ENCFF427WCB (rep 2).

### Result: LOW BUT NON-ZERO accessibility across the element

Viewed at chr11:66,256,955-66,257,171 (element exactly). All six
tracks show signal rising from the left edge and plateauing by
roughly the midpoint. Not baseline. Y-axis maxima are low (one
track reads 8).

**Correction to first impression:** at the wider 1.2 kb view the
element appeared flat, but that was axis compression by the promoter
peak, not absence of signal. Re-checked at element-exact coordinates.

Representative DHS site and cCRE annotation both begin at roughly the
element's right third and extend rightward — so the right end of the
element does fall inside an annotated accessible region; the left
two-thirds do not.

### Interpretation
The element is weakly accessible in motor neurons, far below the
adjacent KLC2 promoter. It is not in closed chromatin, but it is not
a discrete accessible peak either. Consistent with a sequence sitting
at the edge of a strongly accessible promoter rather than functioning
as an independent regulatory element in this cell type.

Does NOT support a strong claim either way about CTCF occupancy here.

### Caveats
- ALS donor-derived only; disease state may alter accessibility
- Signal quantification is visual, not measured. To make any numeric
  claim, download a bigWig and compute mean signal over the element
  vs flanks, as done for H3K27ac.
- Reviewer point 4 (wrong tissue) is now addressed with the correct
  tissue.

## 2026-09-16 — Motor neuron ATAC: full quantification, n=10 experiments / 6 donors

### Supersedes the n=2 result in the 2026-09-15 entry

The previous quantification used three files. ENCODE metadata (scripts/
mn_atac_metadata.py) showed all three came from ONE experiment (ENCSR410DWV),
one biosample (ENCBS435KUD), one donor (ENCDO022IAQ). Worse, ENCFF345PTN is
the pooled file (bio_reps 1,2) and ENCFF576QGF is replicate 1 alone — the
second file is NESTED INSIDE the first, not independent. That is why they
agreed to within 5-7% in every window. The earlier "n=2" was n=1 donor with
one file counted twice, i.e. the same error already flagged for the cortical
interneuron CTCF peaks (GSE117508). Any sentence resting on n=2 is withdrawn.

ENCFF638TCZ was output_type "signal p-value" while the other two were "fold
change over control". The post hoc exclusion is therefore resolved at file
selection under prespecified criterion 2, and no longer requires defending
as a judgment call.

### Additional prespecified rules (set 2026-09-15/16, before full quantification)

6. Where an experiment provides both a pooled file and its constituent
   replicate files, use the pooled file only. Replicate-level files are nested
   within the pooled file and are not independent observations.
7. The unit of replication is the DONOR, not the file or the experiment.
   Average within donor before aggregating across donors.

### Data

ENCODE search for ATAC-seq in motor neuron / spinal cord motor neuron returned
10 experiments across 6 donors; 30 passing bigWigs, reduced to 10 by rule 6.

Quantified by REMOTE RANGE READS over HTTPS (pyBigWig with libcurl). No bigWig
was downloaded; each file is ~1 GB and a full download took ~1 h and failed on
timeout. Four windows totalling ~2.2 kb are fetched per file in seconds.
Verified equivalent to local reads: ENCFF345PTN background 0.553 vs 0.55 and
element 1.244 vs 1.24 against the original local quantification.

### Result (n = 10 experiments, 6 donors; 0 excluded by background floor)

Per-donor medians (range across donors):

  element / left flank   3.20  (3.04 - 3.97)
  element / promoter     0.29  (0.22 - 0.34)
  element / background   2.02  (1.29 - 3.17)

Absolute element signal is the most stable quantity in the table: 1.01-1.38
across all ten files (+/-13%), while background varies 4.3-fold (0.26-1.12)
and promoter 1.5-fold (3.13-4.85). Six independent iPSC lines agree on the
element's absolute accessibility.

### element/background is unstable and is demoted in reporting order

The spread in element/background is driven almost entirely by the denominator.
ENCSR709QRD: background 1.121, element 1.112, el_bkg 0.992 — on that file
alone the element shows NO accessibility above background, while its
element/left-flank is 3.468, in line with every other donor. The background
window (chr11:66,250,000-66,251,000) does not behave as a stable baseline
across files.

Reporting order is therefore element/left-flank first, then element/promoter,
then element/background with the instability noted. All three are reported as
prespecified; only the emphasis is set on the evidence. ENCSR709QRD is retained
and worth citing as a demonstration of why a single dataset is not a result —
the same argument already made against the 2015 nine-cell-line analysis and
against GSE117508.

### Outstanding: window definitions

Background and element windows match the original mn_atac_quant.py exactly.
Left flank and promoter do NOT (left 0.418 vs 0.59; promoter 3.133 vs 2.83 on
ENCFF345PTN) — the rewritten script guessed 500 bp flanking each side. Read the
real coordinates off mn_atac_quant.py, pick ONE definition, rerun, and record
the coordinates here. The el_prom difference (0.29 vs 0.45) is entirely this.
No claim may be written until this is settled.

### New observation, not yet a result

Left flank (~0.35) reads BELOW distant background (~0.55) in most files. The
profile is therefore not a simple ramp: local trough immediately upstream,
element at ~1.15, promoter at ~4.0. Suggestive of a discrete feature but not
decisive.

### Decisive test — still pending

scripts/mn_atac_profile.py: 10-bp-binned profile, chr11:66,255,500-66,258,500,
all ten files overlaid, plus a per-file max-normalised panel so donors of
different depth are comparable. Question unchanged: is there a local minimum
between the element and the promoter summit, reproducible across donors?

Prespecified interpretation is unchanged from 2026-09-15, except that the
prior expectation is revised. On 09-15 I recorded an expectation of monotonic
decay (outcome 2). With six donors agreeing on the element's absolute signal
to +/-13% and a left-flank trough now visible, that expectation is withdrawn
and no outcome is favoured going in.

Peak-call intersection (IDR narrowPeak) still to run as the independent check.

### Standing caveats (unchanged)

- All experiments derive from ALS-patient iPSC lines. KLC2 is not an ALS gene,
  but these are not healthy control motor neurons, and KIF5A — same kinesin-1
  complex — is an ALS gene.
- No DNase-seq exists for this biosample.
- Windows are asymmetric: the promoter abuts the element, so no right flank
  comparable to the left flank can be defined.

### New workstream logged: AlphaGenome (DeepMind), not yet started

Sequence-to-function model, 1 Mb input, single-base-pair resolution across
expression, accessibility, histone marks, TF binding and contact maps. Free
API for non-commercial use. Avsec et al., Nature, Jan 2026,
doi:10.1038/s41586-025-10014-0.

Directly relevant because the project's question is whether a feature inside
the 216 bp can be identified from BASE-RESOLUTION data where bulk chromatin
states cannot resolve it. AlphaGenome is not a reference dataset, so it does
not bear on the claim that no reference dataset resolves the element; it is a
model that generates a testable prediction.

Planned, in order:
  A. Predict the full 216-bp deletion's effect on KLC2 expression across
     tissues. Melo measured +48-74% in fibroblasts and iPSC-MNs, no change in
     whole blood. Does the model reproduce direction and tissue pattern? One
     API call. Decisive about whether to continue. Record the result either way.
  B. In-silico deletion scanning across the 216 bp (20-bp tiles, and the three
     FIMO clusters at 24-52, 81-112, 147-171) to localise the predicted
     repressive effect. This is Aim 1 in silico and would set construct
     boundaries on a prediction rather than on motif rank alone.
  C. Predict the 4-5 bp CTCF core mutation, on expression and predicted CTCF
     binding tracks.

Limitations to state wherever this appears:
  - Prediction, not measurement. Hypothesis-generating; does not replace the
    reporter assay, it tells the reporter assay where to look.
  - AlphaGenome docs state predictions were evaluated on sequences differing
    from reference by relatively small amounts (SNPs, indels) and that large
    differences such as structural variants may be less reliable. A 216-bp
    deletion is at the upper edge of "indel". Indel stitching is supported and
    on by default, so it runs, but the size caveat must be stated. The
    small edits in (B) and (C) are on safer ground than the full deletion in (A).
  - Confirm which tissue ontology terms are available before starting; tibial
    nerve and spinal cord (GTEx) expected, motor neuron uncertain.

Scoping: do (A) only until the poster is assembled. (B) and (C) belong to the
proposal, possibly a Fig 6, and must not delay Fig 5, the framing rewrite, the
Fig 2/3 captions, the methods diagram, or poster assembly.

## 2026-09-16 — AlphaGenome calibration test (Part A): result and reading

### Prespecified before running (recorded at the time, reproduced here)

Melo et al. measured KLC2 +48-74% in patient fibroblasts and iPSC-MNs, and no
change in whole blood. The prediction was to be judged on:
  1. DIRECTION — ALT > REF for KLC2 in fibroblast and nerve/spinal cord.
     Primary criterion.
  2. TISSUE PATTERN — effect in whole blood smaller than in fibroblast/nerve.
  3. MAGNITUDE recorded but not pass/fail; sequence-to-function models are
     more reliable on direction than effect size.
If direction wrong -> report failed calibration, do not proceed to Part B.
If direction right but pattern flat -> proceed with tissue-specificity dropped.

### Method

AlphaGenome API (Avsec et al., Nature 2026, doi:10.1038/s41586-025-10014-0),
alphagenome SDK in its own conda env (`alphagenome`; envs are now klc2 / meme /
alphagenome — recorded because "which env" is what breaks a reanalysis later).
score_variant, RECOMMENDED_VARIANT_SCORERS['RNA_SEQ'], 1-Mb context window
centred on the variant. Deletion encoded VCF-style: anchor base at
chr11:66,256,954, reference = anchor + element, alternate = anchor alone.
Reference bases taken from raw/klc2_del_flank500_hg38.fa, never typed from the
paper — the locus is inverted between builds and hand-entered bases would be
reverse-complemented.

Both interpretations of the deletion were run (del217 as published; del216
dropping the 3' base) because the length discrepancy with Uirá is unresolved.

Output: 57,134 rows x 24 cols, 371 tracks per gene, 77 genes in the window.
processed/ag_klc2_scores.tsv

Aggregation rule: scores averaged per biosample before interpretation. Single
tracks are unreliable — cerebellum gave three tracks spanning -0.010 to +0.012
(q_mean 0.33). The same caution applies to every one-track biosample below,
INCLUDING motor neuron.

### CRITICAL: do not average across lineages

The first summary ranked genes by mean score across all 371 tracks and returned
KLC2 at -0.0125, apparently contradicting the per-tissue result. That number is
meaningless: fibroblast (+0.03) and T cell (-0.33) tracks cancel. Any
cross-lineage mean at this locus describes no tissue. All interpretation below
is within-lineage.

### Result — criterion 1 PASSED

Mean KLC2 raw_score by lineage (del217; del216 within ~15%, see below):

  fibroblast (17 biosamples)        +0.024    quantile ~0.9995
  tibial nerve (2 tracks)           +0.017    quantile  0.9989
  spinal cord (2 tracks)            +0.019    quantile  0.9993
  C1 cervical spinal cord           +0.014    quantile  0.9984
  sciatic nerve                     +0.012    quantile  0.9985
  motor neuron (1 track)            +0.002    quantile  0.8970
  cerebellum (3 tracks)             +0.003    q_mean    0.3326  [uninterpretable]
  lymphoid (20 biosamples)          -0.190    quantile -1.0000

Direction is correct in fibroblast and in every nerve/cord tissue. Raw scores
look small but the quantile scores are 0.9995+ — top hundredth of a percent of
genome-wide variant effects. REPORT THE QUANTILE, NOT THE RAW SCORE.

### Criterion 2 NOT MET — and the failure is the interesting part

Melo measured NO CHANGE in whole blood. The model predicts a strong DECREASE
in lymphoid cells (-0.14 venous blood, -0.19 PBMC, -0.23 NK, -0.28 to -0.33
across all T-cell subsets, quantile -1.0000) — 5-10x larger in magnitude than
the fibroblast increase.

This is not "no change" and must not be written up as confirmation. Criterion 2
as prespecified is NOT MET. Recorded as an unpredicted lineage asymmetry.

It is also a different SHAPE of effect, not merely a sign flip:

  lineage      lncRNA   KLC2    ratio   genes with |mean| > 0.01
  fibroblast   +0.135  +0.024    5.6x   2
  nerve/cord   +0.140  +0.014    9.8x   2
  lymphoid     -0.204  -0.190    1.1x   6  (+ ENSG00000254461, KLC2-AS1,
                                            KLC2-AS2, CNIH2)

Fibroblast and nerve: tight, local, two-gene effect. Lymphoid: the whole
neighbourhood moves together. Worth stating explicitly rather than describing
blood as "the same effect reversed."

### Motor neuron is the weakest positive in the entire set

+0.002, quantile 0.897, ONE track. The cell type the disease is about shows an
order of magnitude less predicted effect than fibroblast. Stated plainly in any
write-up. Single track, so also subject to the aggregation caveat above.

### Specificity test PASSED

Within-lineage gene ranking, 77 genes in the 1-Mb window:

  KLC2 ranks 2 of 77 in fibroblast, nerve/cord AND lymphoid.
  ENSG00000255320 ranks 1 of 77 in all three.
  Everything below the top two drops under 0.010 in fibroblast and nerve.

The deletion moves KLC2 and its immediate neighbours and essentially nothing
else. The prediction is about this gene, not about the locus generally.

### The lncRNA result — open, with a control pending

ENSG00000255320 (the unnamed 13-kb antisense lncRNA whose intron contains the
element) is the TOP-ranked gene: +0.135 fibroblast, +0.140 nerve, 5-10x KLC2.

Note the direction: in fibroblast and nerve the lncRNA and KLC2 move the SAME
way, both up. This is NOT sense/antisense rebalancing (that would show opposite
signs). Both outputs of a shared bidirectional promoter rising when the element
is removed is what removing a repressor acting on that promoter would look like
— consistent with Melo's reporter data and with the current framing, not a
complication of it.

CONFOUND, unresolved: the element lies INSIDE an intron of ENSG00000255320. The
RNA-seq scorer aggregates predicted coverage over the gene body, so deleting
217 bp of the gene's own sequence may move its score for reasons unrelated to
regulation. Until the control below is run, the +0.135 and the rank-1 position
are NOT quotable.

Control design, criteria fixed BEFORE choosing intervals: three 217-bp windows,
fully intronic in ENSG00000255320, >=1 kb from the element, no ENCODE cCRE
overlap, outside the CpG island (66,257,440-66,258,782) and its shore, GC within
~5 points of the element's 54.4%. scripts/fetch_hg38.py + scripts/ag_controls.py.
  controls ~0            -> effect specific to the element; second finding
  controls ~+0.135       -> artifact of deleting within the gene body; drop the
                            lncRNA claim to a footnote, KLC2 result stands alone
  intermediate           -> report the ratio, claim nothing further

THIS DOES NOT REVIVE THE ANTISENSE-LNCRNA MODEL. That model was demoted on GTEx
measurement. A prediction disagreeing with a measurement is something to report,
not something to resolve by preferring the newer result. Write-up status changes
from "killed" to "demoted on GTEx evidence, with a contrary in-silico prediction
noted." Also worth re-reading what the GTEx evidence actually showed — "the
lncRNA is not expressed where the disease is" and "the deletion changes the
lncRNA" are different claims and can both be true.

### 216 vs 217 bp — no longer a concern for this analysis

del216 and del217 agree closely across every tissue (fibroblast +0.033 vs
+0.036 in BJ; lymphoid -0.33 vs -0.32 in Treg; same rank order throughout). The
unresolved base does not change any conclusion here. It still matters for
construct boundaries and for what Uirá submits to ClinVar.

### Coordinate registration — correction to the 2026-09-15 entry

The motif-position check recorded earlier (CTCF motif at index 580, element at
index 500, therefore registration confirmed) is CIRCULAR. The FASTA and the
motif coordinate derive from the same source, so an absolute one-base shift
moves both together and the check still passes. It confirms internal
consistency only, not absolute registration.

Unresolved separately: the hg38 element reads A...C. Given hg19 G at both
interval ends and the build inversion, both hg38 ends should read C. The 3' end
matches; the 5' end does not. Three possibilities: (a) FASTA shifted one base
(0-based start used as 1-based), (b) the symmetric-G property does not survive
the inversion as recorded in the handoff, (c) the hg19 end bases are not both G
and the 216/217 ambiguity rests on something else.

Resolution: scripts/fetch_hg38.py re-fetches the element from UCSC with the same
flanks and diffs against raw/klc2_del_flank500_hg38.fa. Identical => registration
settled. Different => the off-by-one is located. Not a blocker — a +/-1 bp shift
cannot change a prediction made on a 1-Mb window — but it must be settled before
construct design and before Uirá's ClinVar correction.

### Standing limitations for anywhere this appears

- Prediction, not measurement. Hypothesis-generating. Does not replace the
  reporter assay; it tells the reporter assay where to look.
- AlphaGenome docs: predictions evaluated on sequences differing from reference
  by relatively small amounts (SNPs, indels); large differences such as
  structural variants may be less reliable. A 216-bp deletion sits at the upper
  edge of "indel." Indel stitching is supported and on by default.
- Criterion 2 not met (above). Report the discrepancy, do not smooth it.
- Motor neuron near-null, single track.
- Cross-lineage averaging is invalid at this locus.

### Sentence this currently supports (pending the control)

> In-silico prediction (AlphaGenome; Avsec et al. 2026) recapitulates the
> direction of the measured effect: the 216-bp deletion is predicted to increase
> KLC2 in fibroblasts (mean +0.024, 17 biosamples) and in tibial nerve and spinal
> cord, with KLC2 ranked second of 77 genes in the 1-Mb window. The model
> additionally predicts a decrease in lymphoid lineages that was not observed in
> the 2015 whole-blood measurement, and predicts near-null effect in motor
> neurons; both are reported as discrepancies.

### Scoping decision

IN SCOPE NOW (~1 h): element refetch/registration check; three control
deletions; this entry.

DEFERRED until the poster is assembled: in-silico tiling across the 216 bp
(Part B), CTCF core mutation prediction (Part C), any modality beyond RNA_SEQ
(ATAC, CHIP_TF), anything further on the lncRNA beyond the one control.

REASON: still open from the revision plan — Fig 5 rebuild, framing rewrite
through poster/proposal/abstract, Fig 2 and Fig 3 captions, methods flow
diagram, poster assembly, send revised figures to Uirá. The motor neuron ATAC
profile plot (scripts/mn_atac_profile.py) is also still unrun, and unlike
AlphaGenome it decides an actual claim about a measurement. Three sessions have
now gone to a workstream that did not exist yesterday and none of those seven
items has moved.

PLACEMENT: one or two sentences in the proposal's preliminary data plus a
limitations line. NOT a poster figure — a prediction beside six
measurement-based figures invites "why do you trust it," and that answer takes
longer than a poster conversation allows. Have it ready if asked.

## 2026-09-17 — MAJOR CORRECTION: element coordinates were 131 bp off

### Finding

Published hg19 interval chr11:66,024,557-66,024,773 (217 bp, ends G/G) matches
hg38 chr11:66,257,086-66,257,302 at 217/217, FORWARD strand.

The interval used throughout this project, chr11:66,256,955-66,257,171, is
131 bp upstream. Overlap with the true element: 86 bp of 217 (40%).

Verified by direct sequence retrieval from the UCSC REST API (getData/sequence)
for both builds, independently of liftover and of BLAT.
scripts/check_hg19_hg38.py

### Two prior conclusions overturned

1. THE LOCUS IS NOT INVERTED BETWEEN BUILDS. A forward match excludes it. Every
   statement of "hg19 plus strand = hg38 minus strand" in the handoff and in the
   figure captions is wrong.

2. THE LIFTOVER WAS CORRECT. The handoff records a liftover result rejected for
   being "131 bp off" in favour of a BLAT verification. That is exactly this
   offset. The liftover was right; the BLAT result was accepted over it in
   error. Whatever query sequence produced the BLAT hit at 66,256,955 was not
   the sequence at the published hg19 coordinates. To be resolved with Uirá:
   does the 2015 supplementary sequence match the 2015 published coordinates?

This also resolves the A...C anomaly logged 2026-09-16. The true interval reads
G/G in both builds, consistent with the published G at both ends and with the
216/217 ambiguity as described. The old interval read A...C because it began
131 bp early.

### Immediate consequence: the CTCF claim

The CTCF motif at chr11:66,257,035-66,257,065 ends 21 bp BEFORE the true element
begins. It lies entirely outside the deletion.

This removes the basis for Fig 5, the proposal title, Aim 2, and the CTCF
ChIP-qPCR design. Nothing rewritten until FIMO is rerun on the correct sequence.

### Affected, to redo at chr11:66,257,086-66,257,302

  FIMO motif scan (CTCF rank, 14/30 zinc fingers, clusters 24-52/81-112/147-171)
  Fig 5 and the entire CTCF line of argument
  Motor neuron ATAC quantification (element window; promoter and background
    windows were defined relative to the old element and need re-checking)
  AlphaGenome runs (2026-09-16 entry) - wrong deletion scored throughout
  Deletion junction reconstruction and scan
  Locus annotation: TSS distances, cCRE gap, CpG island shore (island starts
    66,257,440 - now 138 bp downstream, not 268), lncRNA intron position,
    Multiz conservation
  ChromHMM (Fig 3) and histone signal (Fig 4) - windows shift 131 bp; the
    resolution conclusion likely survives but must be re-derived, not assumed
  Figs 1-5: highlight position
  REST, cohesin, Hi-C, DNase negatives

### Unaffected

  GTEx reinterpretation (blood 2.8 TPM vs 19-24 in disease tissues) - no
    coordinates involved
  The resolution argument - a 217-bp element ~100 bp from an active TSS is not
    resolvable by bulk chromatin data regardless of a 131-bp shift
  Uirá's framing correction: the element is an established repressive sequence
  All prespecification and reporting discipline

### Old files retained

raw/klc2_del_flank500_hg38.fa and all processed/ outputs derived from it are
KEPT, not overwritten. New work goes to *_TRUE_* filenames. The error and its
detection stay in the record.

### How it was caught

The hg38 element read A...C where the published hg19 ends (G/G) implied C...C
under the assumed inversion. No single-base shift produced C...C, which meant
the problem was not an extraction offset. Checking the hg19 sequence directly
rather than the hg38 registration found the 131-bp offset. A two-base anomaly
was the only visible symptom.

### Next

1. Rebuild element FASTA at true coordinates
2. Rerun FIMO - determines whether any motif story exists and Aim 2's fate
3. Recompute locus annotation
4. Re-quantify ATAC, rerun AlphaGenome, rebuild junction
5. Email Uirá today: he is checking chromatograms for a ONE-base discrepancy and
   needs to know it may be a 131-bp one. Bears on the ClinVar/OMIM corrections
   he offered to make.
## 2026-09-17 (cont.) — Motif analysis at corrected coordinates: three negatives

### FIMO on the corrected element

chr11:66,257,086-66,257,302. FIMO (MEME 5.3.0), JASPAR2024 CORE vertebrates
non-redundant, p < 1e-4. 383 hits in the flank file; 73 inside the 217 bp.

Element GC is 57.6%, NOT the 54.4% recorded in the handoff — that figure was
computed on the wrong sequence. Corrected everywhere.

### Structure found

Three spatial clusters with gaps at 41-74 and 90-102:
  I    16-40    ZNF701, ZNF528, ZNF175, ETS core, ZBTB11
  II   103-141  ETS core, Zbtb2/ZBTB11/ZBTB6, ZNF454, nuclear receptors
  III  149-215  ETS core, NFKB1/2, ZNF768/740/454/460/331/213, INSM1

These SUPERSEDE the old 24-52 / 81-112 / 147-171 clusters, which were on the
wrong sequence.

Two copies of an identical 9-bp ETS core (ACCGGATGT) at element 109-117 and
168-176, both plus strand, 59 bp apart, in unrelated flanks.

CORRECTION LOGGED: I first read this as a 59-bp tandem repeat carrying a
duplicated ETS module. Checked directly — the two 18-mers share only their
first 10 bases (ACCGGATGTG then diverge completely). It is homotypic site
clustering, not a duplication. The reporter design implication differs: mutate
the two cores, not delete two copies of a module.

CTCF still present (p = 2.1e-5) but ~24th rank, inside cluster III. Not
pursued — outside the old element's basis, and Uirá advised dropping it.

### Enrichment test — the analysis Uirá asked for

DESIGN. Element + 30 GC-matched 217-bp background windows in ONE FASTA, ONE
FIMO run, so settings are identical by construction. Background: 57.6% +/- 2.0
GC, non-overlapping, excluding element +/-500 bp and CpG island +/-500 bp.
Sampling regions widened to chr11:66,240,000-66,256,500 and 66,259,500-66,270,000
because at +/-3.0 GC only 14 qualifying windows existed in the original 11 kb —
the element is GC-rich for its neighbourhood, worth recording in itself.
Background windows therefore extend beyond the lncRNA; slight weakening of the
"same intron" control, noted.

Family classification in a SHARED module (scripts/tf_families.py) applied
identically to element and background. Repressor-associated families
(ZNF_C2H2, ZBTB, INSM) PRESPECIFIED before the comparison ran. NFKB excluded
as context-dependent (p50/p52 homodimers repress, heterodimers activate). ETS
excluded as predominantly activating.

Hits collapsed by family+position — JASPAR lists many near-identical matrices
per family and they are not independent observations.

RESULT — THREE NEGATIVES:

  repressor-family hits   element 24 | bg median 19.5 (11-30)
                          8 of 30 bg windows match or exceed | p = 0.290
  total hits              element 45 | bg median 43.5 (21-67)  indistinguishable
  distinct motifs         element 51 | bg median 43.5 (22-82)  p = 0.742
  distinct families       element 8  | bg median 5.0  (4-8)    p = 1.000

>> By motif content, the corrected 217-bp element is indistinguishable from
>> GC-matched sequence in the same region.

This RETIRES the old "14 of 30 FIMO hits are zinc fingers" observation, which
was never tested against a baseline and was in any case computed on the wrong
sequence.

### An intermediate claim raised and killed the same hour

The element showed 0 hits in the "other" family bucket while every background
window had 7-29. I read this as an unusually narrow motif repertoire and said
so. Tested it directly (scripts/family_diversity.py, distinct motifs and
distinct families per window, the former classifier-independent): the element
has MORE families than any background window (8, top of the 4-8 range) and an
ordinary number of distinct motifs. The zero was a classifier artifact — the
eight named families happened to cover everything in that particular sequence.
Claim withdrawn.

Note for figures: the family+position-collapsed table and the raw-hit table
give different counts for the same sequence (element 45 vs 73). Both correct
for their purpose. Do not mix them in one figure.

### Why this strengthens the case

Five independent computational approaches now converge on the same conclusion:
ChromHMM, continuous histone signal, DNase, ATAC, and motif content. The
element is not resolvable from its context by reference data or by sequence
scanning. That convergence is a better argument for the reporter assay than any
single negative.

### STOPPING POINT

The motif line of inquiry is CLOSED. Three negatives from three angles is
sufficient. A fourth variation would be fishing.

### AlphaGenome cluster decomposition and ETS mutations

Run at corrected coordinates. Fibroblast mean KLC2 score:

  full element   217 bp  +0.0339
  cluster I       25 bp  -0.0012
  cluster II      39 bp  +0.0210
  cluster III     67 bp  +0.0178
  ETS core 1       9 bp  +0.0152
  ETS core 2       9 bp  +0.0165

full_217 reproduces the calibration result at corrected coordinates, so that
finding survives the coordinate fix. Cluster I contributes nothing. Clusters II
and III each carry roughly half, and each contains one ETS core. Each 9-bp core
alone reproduces ~45% of the full deletion's predicted effect.

POINT MUTATIONS (ACCGGATGT -> ACCTTATGT; preserves length and spacing, so
isolates the binding site from the geometry):

  core 1 deletion  +0.0152    core 1 mutation  +0.0177
  core 2 deletion  +0.0165    core 2 mutation  +0.0054

Core 1 behaves like a binding site — the mutation reproduces and slightly
exceeds the deletion. Core 2 does NOT — the mutation recovers only a third of
its deletion's effect, so most of what deleting core 2 does comes from removing
sequence, not from destroying the ETS site. Two identical 9-mers behaving
differently argues against a redundant-sites model and implies context matters.

Nerve/cord shows nothing: every edit between -0.008 and +0.004, and full_217 is
only +0.0021 there. The predicted effect is fibroblast-specific. Workable since
fibroblasts are now the primary model, but must be stated, not glossed.

### Controls required before any of the ETS result is quotable

1. Scrambled control at both core positions — show an arbitrary 3-bp change
   does NOT produce the same effect. Without it, "the ETS site matters" is one
   mutation away from "any change there matters."
2. FIMO rescan of the mutant sequences — ACCTTATGT may CREATE a site rather
   than only destroy one. Same trap avoided on the deletion junction.
3. Double mutant — requires predict_sequence on a custom sequence, not
   score_variant. Deferred.

### Scoping

AlphaGenome stops here until the poster/presentation work is done. 12 days to
the 29 Sept meeting and still untouched: UCSC annotation checklist, ATAC
re-quantification at corrected coordinates, ATAC profile plot (written, never
run, and it decides an actual claim about a measurement), junction rebuild, six
figures, the presentation itself.

Next: UCSC checklist (browser work, no code), then ATAC, then figures.

### Uirá's reply received today

Endorsed the GTEx blood reading. Confirmed fibroblasts as a valid primary model
at 19 TPM with patient lines available. Said to stop chasing CTCF — the
architectural hypothesis is closed — and instead scan the 216 bp for repressor
motifs without a prior, then look for footprints in fibroblast ATAC.

The motif scan is done and returned a negative (above). Fibroblast ATAC
footprinting not started — needs BAM files rather than bigWigs, substantially
more work than anything attempted so far. Scope separately, after the meeting.

Meeting offered: Tue 29 Sept or Thu 1 Oct, 17:00 or 18:00 Berlin (11:00/12:00
New York). He asked for a short presentation on background, what has been
tested, what is open, and the hypotheses — and to cover Lucid, his company
(lucid-genomics.com). Read up on Lucid before the call.

Reply drafted, including the 131-bp coordinate finding, since it changes what
he is looking for in the 2015 chromatograms and what the ClinVar/OMIM
corrections should say. Key question to put to him: does the 2015 supplementary
sequence match the 2015 published coordinates?

## 2026-09-17 (cont.) — ATAC at corrected coordinates: a reproducible subpeak
## summit inside the element

### Re-quantification at chr11:66,257,086-66,257,302

Windows redefined relative to the corrected element. NOTE: the OLD promoter
window overlapped the corrected element, so old and new numbers are NOT
comparable and no before/after should be presented.

  background  66,250,000-66,251,000   (unchanged, distant)
  left        66,256,586-66,257,085   (500 bp, immediately upstream)
  element     66,257,086-66,257,302   (217 bp)
  promoter    66,257,303-66,257,802   (500 bp, spans the CpG island start)

Per-donor medians (n = 10 experiments, 6 donors, 0 excluded by the floor):

  element/left        3.11  (2.87-3.39)   <- tightest, and the headline number
  element/promoter    0.29  (0.25-0.34)
  element/background  3.09  (1.76-5.12)   <- still denominator-driven, demoted

Absolute element signal 1.18-2.33 (+/-33%), up from 1.01-1.38 (+/-13%) on the
old interval. The corrected window sits closer to the promoter, so the rise is
expected; the reduced stability is consistent with the window overlapping a
gradient rather than sitting on a flat feature. Observation, not a conclusion.

### Profile test — local maximum over the element's 3' end

scripts/mn_atac_profile.py, chr11:66,255,600-66,258,600, 10-bp bins, all 10
files, plus a per-file max-normalised panel.

A local maximum appears over the element's 3' end, followed by a dip at
~66,257,320-66,257,430, before signal rises into the promoter peak.

QUANTIFIED (scripts/mn_atac_dip.py) rather than eyeballed — and the eye was
wrong. I first read the plot as "reproducible across all ten traces." It is not:

  peak/dip ratio  median 1.17  range 0.87-1.66   >1.0 in 8 of 10

Four files show a clear dip (>=1.3), four are marginal (1.03-1.20), two show
none (0.87, 0.93). The two negatives each have a same-donor replicate well above
1.0, so the split is WITHIN donors, not between them — consistent with
experiment-level variation (depth, S/N) rather than biological absence. That is
an interpretation, not a result.

Window sensitivity (+/-20 bp):
  (66257180,66257300) vs (66257300,66257420)  median 1.19  8/10
  (66257200,66257320) vs (66257320,66257430)  median 1.17  8/10   [primary]
  (66257220,66257340) vs (66257340,66257450)  median 1.12  7/10
Some sensitivity, no cliff. Report the middle window as primary with the range
as the sensitivity check. NOTE the windows were chosen by eye off the plot —
a real weakness, and the reason the peak-call test below matters more.

### Peak calls — the decisive test

TWO SELECTION BUGS CAUGHT AND FIXED, both mine:

1. First pass intersected ~50 peak files for 10 experiments. ENCODE publishes
   conservative IDR, optimal IDR, replicated, pseudoreplicated and IDR-ranked
   files per experiment; the raw intersect counted them as independent. Same
   nesting error as the bigWigs.
2. The download filter (awk on /IDR/) and the selection priority disagreed, so
   for 4 experiments the top-ranked file was never downloaded. Those 4 dropped
   out silently, and the remaining 6 each used whatever ranked highest among
   files that happened to be on disk — mixing peak types across experiments.
   The intermediate "3 of 6" result was therefore not a clean comparison and is
   withdrawn.

FIXED: scripts/fetch_selected_peaks.py picks exactly ONE file per experiment on
a fixed priority (conservative IDR > IDR thresholded > replicated >
pseudoreplicated), downloads it, verifies the gzip, then
scripts/atac_summits.py reads that selection. "IDR ranked peaks" EXCLUDED —
those are unthresholded ranked lists, not called peaks.

RESULT:

  experiments with >=1 summit inside the element:  6 of 10
  donors with >=1:                                 6 of 6
  summit positions: 66,257,254 / 261 / 263 / 267 / 267 / 301

Every donor shows it in at least one experiment. Five of six summits fall within
a 14-bp span. Across a ~1.5 kb accessible region, that clustering is not chance.

IMPORTANT — what this is NOT. Every file calls ONE large accessible region of
~1.5 kb (starts 66,257,013-66,257,047, ends 66,258,550-66,258,880) covering both
the element and the promoter. The element does NOT form a separate peak. What is
reproducible is a SUBPEAK SUMMIT within that domain, inside the element.

The accessible region's 5' boundary sits 40-70 bp upstream of the element start
— the domain begins essentially where the element begins. Worth a figure note.

### CONVERGENCE — three methods, three data types, same ~50 bp

  continuous ATAC profile   local max, 8/10 experiments     element 3' end
  called peak summits       6/10 exp, 6/6 donors            66,257,254-267
  AlphaGenome decomposition ETS core 2 del = +0.017/+0.034  66,257,253-261

ETS core 2 occupies element positions 168-176 = chr11:66,257,253-66,257,261.
Five of six ATAC summits land on it or within 6 bp of it.

This is the first time any dataset has resolved something INSIDE the 217 bp.

### The claim, as it may be written

> Within a single ~1.5 kb accessible region spanning the element and the KLC2
> promoter, a called subpeak summit falls inside the 217-bp element in 6 of 10
> motor neuron ATAC experiments, representing all 6 donors, at
> chr11:66,257,254-66,257,301 (5 of 6 within 14 bp). Continuous signal shows a
> corresponding local maximum in 8 of 10 experiments (median peak/dip 1.17,
> range 0.87-1.66). The summit position coincides with an ETS core at
> 66,257,253-66,257,261 that in-silico deletion independently nominates as
> carrying roughly half the predicted expression effect.

Caveats that travel with it: the element does not form a separate accessible
peak; peak types are not uniform across experiments (record the split from the
"peak types chosen" line); all donors are ALS-derived iPSC lines; AlphaGenome is
prediction, not measurement; the profile windows were chosen by eye.

### Consequence for Aim 1

The construct series is no longer "tile the element." It is: wild type, full
216-bp deletion, ETS core 2 point mutation, ETS core 1 point mutation, both
cores mutated. Four informative conditions, each testing something named, within
the ~$3,000 budget.

This also answers Uirá's question better than the motif scan did. He asked for a
repressor-motif scan without a prior (negative, p = 0.29) and then for
footprints in fibroblast ATAC. The accessibility data he pointed at is what
produced the specific site.

### Priority change

The scrambled-sequence control at both ETS cores has moved from nice-to-have to
LOAD-BEARING. Without it, "the ETS site matters" is one mutation away from "any
change there matters," and the proposal now rests on that site. Same for the
FIMO rescan of the mutant sequences (ACCTTATGT may CREATE a site).

Still ahead of them in order: UCSC annotation checklist, figures, presentation.
