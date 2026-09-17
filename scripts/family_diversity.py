#!/usr/bin/env python3
"""Is the element's motif repertoire narrower than background?
Counts distinct motifs and distinct families per window - classifier-independent
for the motif count, classifier-dependent only for the family count."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tf_families import family
import pandas as pd

d = pd.read_csv("processed/fimo_scanset30/fimo.tsv", sep="\t", comment="#")
d = d.dropna(subset=["motif_alt_id"])
d["family"] = d.motif_alt_id.map(family)

s = (d.groupby("sequence_name")
       .agg(n_hits=("motif_alt_id", "size"),
            n_distinct_motifs=("motif_alt_id", "nunique"),
            n_families=("family", "nunique"))
       .reset_index())

el = s[s.sequence_name == "element"].iloc[0]
bg = s[s.sequence_name != "element"]

print(s.sort_values("n_distinct_motifs").to_string(index=False))
for col in ["n_distinct_motifs", "n_families"]:
    n_le = int((bg[col] <= el[col]).sum())
    print(f"\n{col}: element {el[col]} | background median {bg[col].median():.1f} "
          f"(range {bg[col].min()}-{bg[col].max()})")
    print(f"  background windows at or below element: {n_le} of {len(bg)}"
          f"   empirical p = {(n_le + 1) / (len(bg) + 1):.3f}")
