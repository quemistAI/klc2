#!/usr/bin/env python3
"""Discover all ENCODE motor neuron ATAC-seq experiments and their bigWig files.
Answers: how many distinct donors exist, and which files pass the
prespecified inclusion rule (GRCh38, fold change over control, released)."""

import requests, pandas as pd

H = {"accept": "application/json"}
BASE = "https://www.encodeproject.org"
TERMS = ["motor neuron", "spinal cord motor neuron"]


def get(path, **params):
    params.setdefault("format", "json")
    r = requests.get(f"{BASE}{path}", headers=H, params=params, timeout=60)
    r.raise_for_status()
    return r.json()


# --- find the experiments ---
exp_accs = set()
for term in TERMS:
    try:
        res = get("/search/", type="Experiment", assay_title="ATAC-seq",
                  **{"biosample_ontology.term_name": term},
                  status="released", limit="all")
    except Exception as e:
        print(f"  search failed for '{term}': {e}")
        continue
    hits = res.get("@graph", [])
    print(f"  '{term}': {len(hits)} experiments")
    exp_accs.update(h["accession"] for h in hits)

print(f"\ntotal distinct experiments: {len(exp_accs)}\n")

# --- walk each experiment ---
rows = []
for acc in sorted(exp_accs):
    exp = get(f"/experiments/{acc}/")

    donors, biosamples = set(), set()
    for rep in exp.get("replicates") or []:
        bio = (rep.get("library") or {}).get("biosample") or {}
        if isinstance(bio, str):
            bio = get(bio)
        if bio.get("accession"):
            biosamples.add(bio["accession"])
        dn = bio.get("donor")
        dn = get(dn) if isinstance(dn, str) else (dn or {})
        if dn.get("accession"):
            donors.add(dn["accession"])

    files = get("/search/", type="File", dataset=f"/experiments/{acc}/",
                file_format="bigWig", status="released", limit="all")

    for f in files.get("@graph", []):
        rows.append({
            "experiment": acc,
            "file": f.get("accession"),
            "output_type": f.get("output_type"),
            "assembly": f.get("assembly"),
            "status": f.get("status"),
            "bio_reps": ",".join(map(str, f.get("biological_replicates") or [])),
            "biosamples": ",".join(sorted(biosamples)),
            "donors": ",".join(sorted(donors)),
            "lab": (exp.get("lab") or {}).get("title"),
            "description": (exp.get("description") or "")[:70],
        })

df = pd.DataFrame(rows)
df.to_csv("processed/mn_atac_all_files.tsv", sep="\t", index=False)

# --- apply the prespecified rule ---
passing = df[(df.assembly == "GRCh38") &
             (df.output_type == "fold change over control") &
             (df.status == "released")].copy()
passing.to_csv("processed/mn_atac_passing_files.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
print(passing[["experiment", "file", "bio_reps", "donors"]].to_string(index=False))
print()
print("passing files:      ", len(passing))
print("distinct experiments:", passing.experiment.nunique())
print("distinct donors:     ", len(set(d for s in passing.donors for d in s.split(","))))

