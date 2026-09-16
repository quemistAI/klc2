#!/usr/bin/env python3
"""Resolve ENCODE accessions to direct S3 bigWig URLs (no download)."""

import requests, pandas as pd

H = {"accept": "application/json"}
sel = pd.read_csv("processed/mn_atac_selected.tsv", sep="\t")

rows = []
for _, r in sel.iterrows():
    d = requests.get(f"https://www.encodeproject.org/files/{r.file}/",
                     headers=H, timeout=30).json()
    url = (d.get("cloud_metadata") or {}).get("url") or \
          f"https://www.encodeproject.org{d['href']}"
    rows.append({"experiment": r.experiment, "file": r.file,
                 "donor": r.donors, "url": url})

out = pd.DataFrame(rows)
out.to_csv("processed/mn_atac_urls.tsv", sep="\t", index=False)
print(out.to_string(index=False))
