#!/usr/bin/env python3
"""Find 9-mers that replace the ETS core, preserve GC, destroy the ETS match,
and create no new FIMO hit. Writes candidates for a FIMO screen."""

import itertools, pandas as pd

FA, PAD, N = "raw/element_TRUE_flank500_hg38.fa", 500, 217
CORE1, CORE2 = (109, 117), (168, 176)
WT = "ACCGGATGT"

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
el = seq[PAD:PAD + N]

def gc(s):
    return s.count("G") + s.count("C")

target_gc = gc(WT)
cands = []
for cand in ("".join(p) for p in itertools.product("ACGT", repeat=9)):
    if gc(cand) != target_gc:
        continue
    if "GGA" in cand or "TCC" in cand:      # ETS core, both strands
        continue
    if sum(a != b for a, b in zip(cand, WT)) < 4:   # ensure a real change
        continue
    cands.append(cand)

print(f"{len(cands)} candidate 9-mers; writing a random 300 for screening")
import random
random.seed(23)
pick = random.sample(cands, min(300, len(cands)))

with open("raw/mutant_screen.fa", "w") as f:
    f.write(f">WT\n{el}\n")
    for i, c in enumerate(pick):
        for lbl, (a, b) in (("c1", CORE1), ("c2", CORE2)):
            s = el[:a-1] + c + el[b:]
            f.write(f">{lbl}_{i:03d}_{c}\n{s}\n")

print("wrote raw/mutant_screen.fa")

