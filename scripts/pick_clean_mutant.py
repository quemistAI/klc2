#!/usr/bin/env python3
"""Rank screened mutants: destroys the ETS stack, creates nothing."""

import pandas as pd

d = pd.read_csv("processed/fimo_mutant_screen/fimo.tsv", sep="\t", comment="#")
d = d.dropna(subset=["motif_alt_id"])
d["key"] = d.motif_alt_id + "@" + d.start.astype(str) + "-" + d.stop.astype(str)
wt = set(d[d.sequence_name == "WT"].key)

rows = []
for name, g in d[d.sequence_name != "WT"].groupby("sequence_name"):
    mut = set(g.key)
    rows.append({"name": name, "core": name.split("_")[0],
                 "ninemer": name.split("_")[-1],
                 "lost": len(wt - mut), "gained": len(mut - wt)})

r = pd.DataFrame(rows).sort_values(["core", "gained", "lost"],
                                   ascending=[True, True, False])
r.to_csv("processed/mutant_screen_ranked.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
for core in ["c1", "c2"]:
    print(f"\n=== {core}: best candidates (gained=0 preferred) ===")
    print(r[r.core == core].head(10).to_string(index=False))

best = {}
for core in ["c1", "c2"]:
    sub = r[(r.core == core) & (r.gained == 0)]
    if len(sub):
        best[core] = sub.iloc[0].ninemer
print("\nbest clean 9-mers:", best or "none with gained=0 — take the minimum")
