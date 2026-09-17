#!/usr/bin/env python3
"""List AlphaGenome ontology terms matching the SPOAN-relevant tissues.
RUN IN: conda env `alphagenome`.  Writes to processed/ (shared with klc2 env)."""

import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env

API_KEY = check_env()

import pandas as pd
from alphagenome.models import dna_client

model = dna_client.create(API_KEY)
meta = model.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)

frames = []
for attr in ["rna_seq", "atac", "dnase", "chip_tf", "chip_histone"]:
    df = getattr(meta, attr, None)
    if df is None or len(df) == 0:
        continue
    d = df.copy()
    d["modality"] = attr
    frames.append(d)

all_meta = pd.concat(frames, ignore_index=True)
os.makedirs("processed", exist_ok=True)
all_meta.to_csv("processed/ag_output_metadata.tsv", sep="\t", index=False)
print("modalities:", all_meta.modality.value_counts().to_dict(), "\n")

PAT = re.compile(r"nerve|spinal|fibroblast|blood|motor neuron|neuron|cerebell|cortex|muscle",
                 re.I)
cols = [c for c in all_meta.columns
        if re.search(r"ontology|biosample|name|term|tissue", c, re.I)]
hit = all_meta[all_meta[cols].astype(str).apply(
    lambda r: bool(PAT.search(" ".join(r))), axis=1)]

pd.set_option("display.width", 220)
pd.set_option("display.max_rows", 400)
print(hit[cols + ["modality"]].drop_duplicates().to_string(index=False))
