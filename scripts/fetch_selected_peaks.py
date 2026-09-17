#!/usr/bin/env python3
"""Download exactly one peak file per experiment, chosen by a fixed priority,
so every experiment is compared on the same peak type where available."""

import os, subprocess, pandas as pd

PRIORITY = ["conservative IDR thresholded peaks",
            "IDR thresholded peaks",
            "replicated peaks",
            "pseudoreplicated peaks"]

urls = pd.read_csv("processed/mn_atac_peak_urls.tsv", sep="\t")
urls = urls[urls.output_type.isin(PRIORITY)].copy()
urls["rank"] = urls.output_type.map(PRIORITY.index)

one = (urls.sort_values(["experiment", "rank", "file"])
           .groupby("experiment", as_index=False).first())
one.to_csv("processed/mn_atac_peak_selected.tsv", sep="\t", index=False)
print(one[["experiment", "file", "output_type"]].to_string(index=False))
print("\npeak types chosen:", one.output_type.value_counts().to_dict())

os.makedirs("raw/mn_atac_peaks", exist_ok=True)
for _, r in one.iterrows():
    out = f"raw/mn_atac_peaks/{r.file}.bed.gz"
    if os.path.exists(out) and subprocess.run(["gzip", "-t", out]).returncode == 0:
        print(f"have {r.file}")
        continue
    print(f"fetching {r.file}")
    subprocess.run(["curl", "-sSL", "-o", out, r.url], check=True)
    if subprocess.run(["gzip", "-t", out]).returncode != 0:
        print(f"  CORRUPT after download: {r.file}")
