#!/usr/bin/env python3
"""In-silico deletion of each motif cluster within the corrected element.
RUN IN: conda env `alphagenome`."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env
API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

FA, PAD = "raw/element_TRUE_flank500_hg38.fa", 500
EL_START = 66257086                      # hg38, 1-based

# label -> (start, end) in element coordinates, 1-based inclusive
TARGETS = {
    "full_217":    (1, 217),
    "clusterI":    (16, 40),
    "clusterII":   (103, 141),
    "clusterIII":  (149, 215),
    "ets_core_1":  (109, 117),
    "ets_core_2":  (168, 176),
}

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

frames = []
for label, (a, b) in TARGETS.items():
    g_start = EL_START + a - 1
    core = seq[PAD + a - 1: PAD + b]
    anchor = seq[PAD + a - 2]
    v = genome.Variant(chromosome="chr11", position=g_start - 1,
                       reference_bases=anchor + core, alternate_bases=anchor)
    print(f"{label}: {len(core)} bp at chr11:{g_start}-{EL_START + b - 1}")
    s = model.score_variant(
        interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=v, variant_scorers=[scorer])
    df = variant_scorers.tidy_scores(s)
    df["deletion"] = label
    frames.append(df)

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_cluster_scores.tsv", sep="\t", index=False)

FIB = r"fibroblast|BJ|IMR-90|HFFc6|AG04450"
NERVE = r"tibial nerve|sciatic nerve|spinal cord|motor neuron"
for name, pat in (("fibroblast", FIB), ("nerve/cord", NERVE)):
    sub = out[out.biosample_name.astype(str).str.contains(pat, case=False, na=False)]
    piv = (sub[sub.gene_name.isin(["KLC2", "ENSG00000255320"])]
           .groupby(["deletion", "gene_name"]).raw_score.mean().unstack().round(4))
    print(f"\n=== {name} mean score ===")
    print(piv.to_string())
