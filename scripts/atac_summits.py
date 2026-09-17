#!/usr/bin/env python3
"""Do called ATAC peak summits fall inside the corrected SPOAN element?
One peak file per experiment; aggregate to donor."""

import glob, gzip, os, pandas as pd

EL = (66257086, 66257302)
one = pd.read_csv("processed/mn_atac_peak_selected.tsv", sep="\t")
sel = pd.read_csv("processed/mn_atac_selected.tsv", sep="\t")[["experiment", "donors"]]
one = one.merge(sel, on="experiment", how="left")
print(one[["experiment", "file", "output_type", "donors"]].to_string(index=False))
print(one[["experiment", "file", "output_type", "donors"]].to_string(index=False))

rows = []
for _, r in one.iterrows():
    path = f"raw/mn_atac_peaks/{r.file}.bed.gz"
    if not os.path.exists(path):
        print(f"  MISSING {r.file}"); continue
    try:
        with gzip.open(path, "rt") as fh:
            summits = []
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] != "chr11":
                    continue
                start, end = int(f[1]), int(f[2])
                if end < EL[0] or start > EL[1]:
                    continue
                s = start + int(f[9])
                summits.append((s, float(f[6])))
    except Exception as e:
        print(f"  FAILED {r.file}: {e}"); continue

    inside = [s for s, _ in summits if EL[0] <= s <= EL[1]]
    rows.append({"experiment": r.experiment, "donor": r.donors, "file": r.file,
                 "n_summits": len(summits), "n_inside": len(inside),
                 "inside_positions": ",".join(map(str, sorted(inside))) or "-"})

d = pd.DataFrame(rows)
d.to_csv("processed/atac_summits_TRUE.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print("\n" + d.to_string(index=False))

print(f"\nexperiments with >=1 summit inside the element: "
      f"{(d.n_inside > 0).sum()} of {len(d)}")
print(f"donors with >=1: {d[d.n_inside > 0].donor.nunique()} of {d.donor.nunique()}")
allpos = [int(p) for s in d.inside_positions for p in s.split(",") if p != "-"]
if allpos:
    print(f"summit positions: {min(allpos)}-{max(allpos)}  "
          f"(element {EL[0]}-{EL[1]}; ETS core 2 at 66257253-66257261)")
