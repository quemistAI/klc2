#!/usr/bin/env python3
"""AlphaGenome calibration at CORRECTED coordinates (chr11:66,257,086-66,257,302).
Reruns the tissue-pattern and gene-specificity results that were previously
computed on the wrong 131-bp-shifted interval.
RUN IN: conda env `alphagenome`."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env
API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

FA, PAD = "raw/element_TRUE_flank500_hg38.fa", 500
EL_START, EL_END = 66257086, 66257302
N = EL_END - EL_START + 1

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
el = seq[PAD:PAD + N]
assert len(el) == 217, f"element is {len(el)} bp"
anchor = seq[PAD - 1]
print(f"element {N} bp, ends {el[0]}/{el[-1]}  (expect G/G)")
print(f"5' {el[:20]}  3' {el[-20:]}")

VARIANTS = {
    "del217": genome.Variant(chromosome="chr11", position=EL_START - 1,
                             reference_bases=anchor + el, alternate_bases=anchor),
    "del216": genome.Variant(chromosome="chr11", position=EL_START - 1,
                             reference_bases=anchor + el[:-1], alternate_bases=anchor),
}

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

frames = []
for label, v in VARIANTS.items():
    print(f"\nscoring {label} ({len(v.reference_bases)-1} bp deleted) ...")
    s = model.score_variant(
        interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=v, variant_scorers=[scorer])
    df = variant_scorers.tidy_scores(s)
    df["variant"] = label
    frames.append(df)

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_calibration_TRUE.tsv", sep="\t", index=False)
print(f"\nwrote processed/ag_calibration_TRUE.tsv  ({out.shape[0]} rows)")
