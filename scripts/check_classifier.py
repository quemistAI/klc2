#!/usr/bin/env python3
"""Which motifs land in 'other', and do element and background differ?"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tf_families import family

import pandas as pd

d = pd.read_csv("processed/fimo_scanset30/fimo.tsv", sep="\t", comment="#")
d = d.dropna(subset=["motif_alt_id"])
d["family"] = d.motif_alt_id.map(family)
d["is_el"] = d.sequence_name == "element"

oth = d[d.family == "other"]
print("=== motifs classed 'other' ===")
print(f"element: {oth.is_el.sum()} hits | background: {(~oth.is_el).sum()} hits\n")
print("top unclassified motif names (background):")
print(oth[~oth.is_el].motif_alt_id.value_counts().head(25).to_string())
print("\ntop unclassified motif names (element):")
print(oth[oth.is_el].motif_alt_id.value_counts().head(25).to_string())

