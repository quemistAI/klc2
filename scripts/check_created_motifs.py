#!/usr/bin/env python3
"""Which motifs does each mutation destroy, and which does it CREATE?"""

import pandas as pd

d = pd.read_csv("processed/fimo_ets_mutants/fimo.tsv", sep="\t", comment="#")
d = d.dropna(subset=["motif_alt_id"])
d["key"] = d.motif_alt_id + "@" + d.start.astype(str) + "-" + d.stop.astype(str)

wt = set(d[d.sequence_name == "WT"].key)
defs = pd.read_csv("processed/ets_mutant_defs.tsv", sep="\t").set_index("name")

for name in d.sequence_name.unique():
    if name == "WT":
        continue
    mut = set(d[d.sequence_name == name].key)
    lost, gained = wt - mut, mut - wt
    r = defs.loc[name] if name in defs.index else None
    print(f"\n=== {name}  ({r.wt} -> {r.mut} at {r.el_start}-{r.el_stop})"
          if r is not None else f"\n=== {name}")
    print(f"  lost {len(lost)}, gained {len(gained)}")
    if gained:
        g = d[(d.sequence_name == name) & (d.key.isin(gained))]
        print("  CREATED:")
        print(g[["motif_alt_id", "start", "stop", "p-value",
                 "matched_sequence"]].sort_values("p-value").head(10)
              .to_string(index=False))
