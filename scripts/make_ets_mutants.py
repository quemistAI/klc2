#!/usr/bin/env python3
"""Build WT + mutant element sequences for FIMO rescan and AlphaGenome scoring.
Writes raw/ets_mutants.fa (217 bp each) and processed/ets_mutant_defs.tsv."""

import pandas as pd

FA, PAD, N = "raw/element_TRUE_flank500_hg38.fa", 500, 217
CORE1, CORE2 = (109, 117), (168, 176)      # element coords, 1-based inclusive
WT_CORE = "ACCGGATGT"

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
el = seq[PAD:PAD + N]
assert len(el) == N

for lbl, (a, b) in (("core1", CORE1), ("core2", CORE2)):
    assert el[a-1:b] == WT_CORE, f"{lbl}: found {el[a-1:b]}"

def sub(s, a, b, new):
    assert len(new) == b - a + 1
    return s[:a-1] + new + s[b:]

def flip3(s):                       # change 3 bases, preserve length
    m = {"A": "C", "C": "A", "G": "T", "T": "G"}
    return "".join(m[c] for c in s[:3]) + s[3:]

defs, recs = [], [("WT", el)]
for lbl, (a, b) in (("core1", CORE1), ("core2", CORE2)):
    nb_a, nb_b = b + 1, b + 9                      # 9 bp immediately 3'
    variants = {
        f"{lbl}_etsmut":   (a, b, "ACCTTATGT"),
        f"{lbl}_scramble": (a, b, WT_CORE[::-1]),  # reversed: composition identical
        f"{lbl}_neighbor": (nb_a, nb_b, flip3(el[nb_a-1:nb_b])),
    }
    for name, (x, y, new) in variants.items():
        recs.append((name, sub(el, x, y, new)))
        defs.append({"name": name, "el_start": x, "el_stop": y,
                     "wt": el[x-1:y], "mut": new})

# both cores at once
both = sub(sub(el, *CORE1, "ACCTTATGT"), *CORE2, "ACCTTATGT")
recs.append(("both_etsmut", both))
defs.append({"name": "both_etsmut", "el_start": CORE1[0], "el_stop": CORE2[1],
             "wt": "two cores", "mut": "ACCTTATGT x2"})

with open("raw/ets_mutants.fa", "w") as f:
    for name, s in recs:
        assert len(s) == N, f"{name} is {len(s)} bp"
        f.write(f">{name}\n")
        for i in range(0, N, 60):
            f.write(s[i:i+60] + "\n")

d = pd.DataFrame(defs)
d.to_csv("processed/ets_mutant_defs.tsv", sep="\t", index=False)
print(d.to_string(index=False))
print(f"\nwrote raw/ets_mutants.fa: {len(recs)} sequences x {N} bp")
