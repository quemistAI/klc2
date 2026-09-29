#!/usr/bin/env python3
"""Score every ETS mutant and control against the WT element.
RUN IN: conda env `alphagenome`."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env
API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

FA, PAD = "raw/element_TRUE_flank500_hg38.fa", 500
EL_START = 66257086

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
defs = pd.read_csv("processed/ets_mutant_defs.tsv", sep="\t")
defs = defs[defs.name != "both_etsmut"]        # needs predict_sequence, deferred

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

frames = []
for _, r in defs.iterrows():
    a, b = int(r.el_start), int(r.el_stop)
    ref = seq[PAD + a - 1: PAD + b]
    assert ref == r.wt, f"{r['name']}: expected {r.wt}, found {ref}"
    v = genome.Variant(chromosome="chr11", position=EL_START + a - 1,
                       reference_bases=ref, alternate_bases=r.mut)
    print(f"scoring {r['name']} ({ref} -> {r.mut}) ...")
    s = model.score_variant(
        interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=v, variant_scorers=[scorer])
    df = variant_scorers.tidy_scores(s)
    df["edit"] = r["name"]
    frames.append(df)

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_ets_control_scores.tsv", sep="\t", index=False)

FIB = r"fibroblast|BJ|IMR-90|HFFc6|AG04450"
sub = out[out.biosample_name.astype(str).str.contains(FIB, case=False, na=False)]
piv = (sub[sub.gene_name.isin(["KLC2", "ENSG00000255320"])]
       .groupby(["edit", "gene_name"]).raw_score.mean().unstack().round(4))
pd.set_option("display.width", 200)
print("\n=== fibroblast mean score ===")
print(piv.to_string())
