#!/usr/bin/env python3
"""Build a FASTA with the corrected element plus 10 GC-matched background
windows from the same lncRNA intron. Stdlib only - runs in any env."""

import json, random, subprocess, sys

random.seed(17)                      # fixed seed = reproducible selection
CHROM, WIDTH = "chr11", 217
EL_START, EL_END = 66257086, 66257302
TARGET_GC, TOL, N_WANTED = 57.6, 2.0, 30

# sample from the lncRNA region, excluding the element +/-500 bp and the
# CpG island (66,257,440-66,258,782) plus a 500 bp margin
REGIONS = [(66240000, 66256500), (66259500, 66270000)]
EXCLUDE = [(EL_START - 500, EL_END + 500), (66256940, 66259282)]


def fetch(start1, end1):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
           f"chrom={CHROM};start={start1-1};end={end1}")
    out = subprocess.run(["curl", "-sSL", url], capture_output=True, text=True).stdout
    return json.loads(out)["dna"].upper()


def gc(s):
    return 100 * (s.count("G") + s.count("C")) / len(s)


def blocked(s, e):
    return any(not (e < a or s > b) for a, b in EXCLUDE)


element = fetch(EL_START, EL_END)
assert len(element) == WIDTH, f"element is {len(element)} bp"
print(f"element: {WIDTH} bp, GC {gc(element):.1f}%")

# pull each region once, then slice locally (avoids hundreds of API calls)
pool = {}
for rs, re_ in REGIONS:
    pool[(rs, re_)] = fetch(rs, re_)
    print(f"fetched {CHROM}:{rs}-{re_}  ({re_-rs+1} bp)")

picked, tries = [], 0
while len(picked) < N_WANTED and tries < 20000:
    tries += 1
    (rs, re_), seq = random.choice(list(pool.items()))
    off = random.randint(0, len(seq) - WIDTH)
    s, e = rs + off, rs + off + WIDTH - 1
    if blocked(s, e):
        continue
    w = seq[off:off + WIDTH]
    if "N" in w or abs(gc(w) - TARGET_GC) > TOL:
        continue
    if any(abs(s - p[0]) < WIDTH for p in picked):    # non-overlapping
        continue
    picked.append((s, e, w))

print(f"search: {tries} attempts for {len(picked)} windows")
if len(picked) < N_WANTED:
    sys.exit(f"only found {len(picked)} matching windows; widen TOL or REGIONS")

with open("raw/scan_set.fa", "w") as f:
    f.write(f">element {CHROM}:{EL_START}-{EL_END} GC={gc(element):.1f}\n{element}\n")
    for i, (s, e, w) in enumerate(sorted(picked), 1):
        f.write(f">bg{i:02d} {CHROM}:{s}-{e} GC={gc(w):.1f}\n{w}\n")

print(f"\nwrote raw/scan_set.fa: 1 element + {len(picked)} background windows")
for i, (s, e, w) in enumerate(sorted(picked), 1):
    print(f"  bg{i:02d}  {CHROM}:{s}-{e}  GC {gc(w):.1f}%")

