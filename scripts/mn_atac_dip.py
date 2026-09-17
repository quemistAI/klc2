#!/usr/bin/env python3
"""Quantify the local maximum over the element and the dip between it and the
promoter peak. Eyeballing a plot is not a result."""

import glob, numpy as np

cand = sorted(glob.glob("processed/mn_atac_profile*.npz"))
if not cand:
    raise SystemExit("no profile .npz found - run scripts/mn_atac_profile.py first")
path = cand[-1]
print(f"using {path}\n")

d = np.load(path)
x = d["x"]

PEAK = (66257200, 66257320)   # apparent local max, 3' end of element
DIP  = (66257320, 66257430)   # apparent minimum
PROM = (66257600, 66257800)   # promoter peak

rows = []
for f in d.files:
    if f == "x":
        continue
    v = d[f]
    pk  = v[(x >= PEAK[0]) & (x < PEAK[1])].max()
    dip = v[(x >= DIP[0])  & (x < DIP[1])].min()
    prm = v[(x >= PROM[0]) & (x < PROM[1])].max()
    rows.append((f, pk, dip, prm, pk / dip if dip else np.nan))
    print(f"{f}  peak {pk:6.2f}  dip {dip:6.2f}  ratio {pk/dip:5.2f}  prom {prm:6.2f}")

r = np.array([x[4] for x in rows])
print(f"\nn files {len(r)}")
print(f"peak/dip ratio  median {np.median(r):.2f}  range {r.min():.2f}-{r.max():.2f}")
print(f"files with ratio > 1.0: {(r > 1.0).sum()} of {len(r)}")
