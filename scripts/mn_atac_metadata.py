#!/usr/bin/env python3
"""Pull ENCODE metadata for motor neuron ATAC-seq files.
Fills in the file-selection table for the prespecified inclusion rule."""

import requests, pandas as pd, sys

H = {"accept": "application/json"}
BASE = "https://www.encodeproject.org"

FILES = [
    "ENCFF345PTN",
    "ENCFF576QGF",
    "ENCFF638TCZ",
    # <- add the other 7 accessions here
]

def get(path):
    r = requests.get(f"{BASE}{path}", headers=H, timeout=30)
    r.raise_for_status()
    return r.json()

rows = []
for f in FILES:
    try:
        d = get(f"/files/{f}/")
    except Exception as e:
        print(f"{f}: FAILED {e}", file=sys.stderr)
        continue

    ds = get(d["dataset"]) if isinstance(d.get("dataset"), str) else d.get("dataset", {})

    donor = biosample = term = None
    reps = ds.get("replicates") or []
    if reps:
        bio = ((reps[0].get("library") or {}).get("biosample")) or {}
        if isinstance(bio, str):
            bio = get(bio)
        biosample = bio.get("accession")
        term = (bio.get("biosample_ontology") or {}).get("term_name")
        dn = bio.get("donor")
        donor = get(dn).get("accession") if isinstance(dn, str) else (dn or {}).get("accession")

    rows.append({
        "file": f,
        "output_type": d.get("output_type"),
        "file_format": d.get("file_format"),
        "assembly": d.get("assembly"),
        "status": d.get("status"),
        "experiment": ds.get("accession"),
        "biosample": biosample,
        "term_name": term,
        "donor": donor,
        "bio_reps": ",".join(map(str, d.get("biological_replicates") or [])),
        "tech_reps": ",".join(map(str, d.get("technical_replicates") or [])),
    })

df = pd.DataFrame(rows)
df.to_csv("processed/mn_atac_metadata.tsv", sep="\t", index=False)
print(df.to_string(index=False))
print()
print("distinct experiments:", df.experiment.nunique())
print("distinct donors:", df.donor.nunique())
print("output types:", df.output_type.value_counts().to_dict())
