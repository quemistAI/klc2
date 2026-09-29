#!/usr/bin/env python3
"""Select matched 217-bp control deletions in ENSG00000255320 introns.

Prespecified criteria (fixed 2026-09-29, before any scoring):
  - exactly 217 bp
  - within the ENSG00000255320 span, overlapping NO exon of ANY isoform
  - >= 1 kb from the SPOAN element
  - no ENCODE cCRE overlap
  - outside the CpG island 66,257,440-66,258,782 plus a 500 bp margin
  - GC within 5 points of the element's 57.6%
Stdlib + requests. Run in klc2."""

import json, subprocess, sys

CHROM = "chr11"
GENE_START, GENE_END = 66244549, 66258438      # widest isoform span
EL_START, EL_END = 66257086, 66257302
CPG = (66257440, 66258782)
WIDTH, TARGET_GC, TOL = 217, 57.6, 5.0
MIN_DIST = 1000
CPG_MARGIN = 500


def api(track, start, end):
    url = (f"https://api.genome.ucsc.edu/getData/track?genome=hg38;"
           f"track={track};chrom={CHROM};start={start};end={end}")
    out = subprocess.run(["curl", "-sSL", url], capture_output=True, text=True).stdout
    try:
        return json.loads(out)
    except Exception:
        sys.exit(f"UCSC returned non-JSON for {track}:\n{out[:300]}")


def seqfetch(start1, end1):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
           f"chrom={CHROM};start={start1-1};end={end1}")
    out = subprocess.run(["curl", "-sSL", url], capture_output=True, text=True).stdout
    return json.loads(out)["dna"].upper()


# --- exons of every overlapping transcript -----------------------------------
d = api("knownGene", GENE_START - 5000, GENE_END + 5000)
items = d.get("knownGene") or next(v for k, v in d.items() if isinstance(v, list))

exons = []
for t in items:
    if t.get("strand") != "-":
        continue
    if t.get("chromEnd", 0) < GENE_START or t.get("chromStart", 1e12) > GENE_END:
        continue
    starts = t.get("chromStarts") or t.get("exonStarts") or ""
    sizes  = t.get("blockSizes")  or ""
    if isinstance(starts, str) and "," in starts:
        base = t.get("chromStart", 0)
        ss = [int(x) for x in starts.strip(",").split(",") if x]
        zz = [int(x) for x in str(sizes).strip(",").split(",") if x]
        for s, z in zip(ss, zz):
            a = base + s if max(ss) < 1e6 else s     # relative vs absolute
            exons.append((a + 1, a + z))
print(f"collected {len(exons)} exon intervals from overlapping minus-strand transcripts")

# --- cCREs -------------------------------------------------------------------
c = api("encodeCcreCombined", GENE_START - 5000, GENE_END + 5000)
citems = c.get("encodeCcreCombined") or next((v for v in c.values() if isinstance(v, list)), [])
ccres = [(x["chromStart"] + 1, x["chromEnd"]) for x in citems]
print(f"collected {len(ccres)} cCREs in the window")

blocked = exons + ccres + [(EL_START - MIN_DIST, EL_END + MIN_DIST),
                           (CPG[0] - CPG_MARGIN, CPG[1] + CPG_MARGIN)]


def clash(s, e):
    return any(not (e < a or s > b) for a, b in blocked)


# --- scan ---------------------------------------------------------------------
region = seqfetch(GENE_START, GENE_END)
cands = []
pos = GENE_START
while pos + WIDTH - 1 <= GENE_END:
    s, e = pos, pos + WIDTH - 1
    if clash(s, e):
        pos += 50
        continue
    w = region[s - GENE_START: e - GENE_START + 1]
    if "N" in w:
        pos += 50
        continue
    gc = 100 * (w.count("G") + w.count("C")) / WIDTH
    if abs(gc - TARGET_GC) <= TOL:
        cands.append((s, e, round(gc, 1)))
    pos += 50

print(f"\n{len(cands)} candidate windows pass all criteria")
if not cands:
    sys.exit("none found - widen TOL or relax MIN_DIST, and record that in the notebook")

# spread the picks out: nearest-GC candidate from each third of the region
picks, seen = [], []
thirds = [(GENE_START, GENE_START + (GENE_END - GENE_START) // 3),
          (GENE_START + (GENE_END - GENE_START) // 3,
           GENE_START + 2 * (GENE_END - GENE_START) // 3),
          (GENE_START + 2 * (GENE_END - GENE_START) // 3, GENE_END)]
for lo, hi in thirds:
    inthird = [c for c in cands if lo <= c[0] < hi
               and all(abs(c[0] - p[0]) > 2 * WIDTH for p in picks)]
    if inthird:
        picks.append(min(inthird, key=lambda c: abs(c[2] - TARGET_GC)))

print("\nselected controls:")
for i, (s, e, gc) in enumerate(picks, 1):
    print(f"  ctrl{i}  {CHROM}:{s:,}-{e:,}  GC {gc}%  "
          f"({abs(s - EL_START):,} bp from element)")

print("\nfetch commands:")
for i, (s, e, gc) in enumerate(picks, 1):
    print(f"python3 scripts/fetch_hg38.py {CHROM} {s} {e} "
          f"raw/lnc_ctrl{i}_flank500_hg38.fa")

print("\npaste into scripts/ag_lnc_controls.py:")
print('DELS = {')
print('    "element": ("raw/element_TRUE_flank500_hg38.fa", 66257086, 66257302),')
for i, (s, e, gc) in enumerate(picks, 1):
    print(f'    "ctrl{i}":   ("raw/lnc_ctrl{i}_flank500_hg38.fa", {s}, {e}),')
print('}')
