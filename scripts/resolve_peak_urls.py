#!/usr/bin/env python3
"""Resolve IDR-thresholded narrowPeak URLs for the motor neuron ATAC experiments."""

import requests, pandas as pd

H = {"accept": "application/json"}
BASE = "https://www.encodeproject.org"
sel = pd.read_csv("processed/mn_atac_selected.tsv", sep="\t")

rows = []
for exp in sel.experiment.unique():
    r = requests.get(f"{BASE}/search/", headers=H, timeout=60,
                     params={"type": "File", "dataset": f"/experiments/{exp}/",
                             "file_format": "bed", "status": "released",
                             "format": "json", "limit": "all"}).json()
    for f in r.get("@graph", []):
        if f.get("assembly") != "GRCh38":
            continue
        if "peaks" not in (f.get("output_type") or ""):
            continue
        rows.append({"experiment": exp, "file": f["accession"],
                     "output_type": f["output_type"],
                     "url": BASE + f["href"]})

out = pd.DataFrame(rows)
out.to_csv("processed/mn_atac_peak_urls.tsv", sep="\t", index=False)
print(out.to_string(index=False))
print("\noutput types found:", out.output_type.value_counts().to_dict())#!/usr/bin/env python3
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
