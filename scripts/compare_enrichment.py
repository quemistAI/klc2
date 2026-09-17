#!/usr/bin/env python3
"""Is the corrected SPOAN element enriched for repressor-family motifs
relative to GC-matched windows from the same intron?"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tf_families import family, REPRESSOR_FAMILIES

import pandas as pd

d = pd.read_csv("processed/fimo_scanset30/fimo.tsv", sep="\t", comment="#")
d = d.dropna(subset=["motif_alt_id"])
d["family"] = d.motif_alt_id.map(family)

# collapse: one count per (sequence, family, position) - JASPAR lists many
# near-identical matrices per family and they are not independent hits
u = d.drop_duplicates(["sequence_name", "family", "start", "stop"])

tab = (u.groupby(["sequence_name", "family"]).size()
         .unstack(fill_value=0))
tab["REPRESSOR_TOTAL"] = tab[[c for c in tab.columns
                              if c in REPRESSOR_FAMILIES]].sum(axis=1)
tab["ALL_HITS"] = tab.drop(columns=["REPRESSOR_TOTAL"]).sum(axis=1)

pd.set_option("display.width", 200)
print(tab.to_string())

el = tab.loc["element"]
bg = tab.drop(index="element")

print("\n=== repressor-family hits ===")
print(f"element:    {el.REPRESSOR_TOTAL}")
print(f"background: median {bg.REPRESSOR_TOTAL.median():.1f}  "
      f"range {bg.REPRESSOR_TOTAL.min()}-{bg.REPRESSOR_TOTAL.max()}  "
      f"mean {bg.REPRESSOR_TOTAL.mean():.1f}")
n_ge = int((bg.REPRESSOR_TOTAL >= el.REPRESSOR_TOTAL).sum())
print(f"background windows matching or exceeding element: {n_ge} of {len(bg)}")
print(f"empirical p = {(n_ge + 1) / (len(bg) + 1):.3f}")

print("\n=== all hits (specificity check) ===")
print(f"element {el.ALL_HITS}   background median {bg.ALL_HITS.median():.1f} "
      f"({bg.ALL_HITS.min()}-{bg.ALL_HITS.max()})")
