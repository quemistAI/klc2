#!/usr/bin/env python3
"""Clean ETS-core replacement test: 9-mers that destroy the ETS stack and
create no new FIMO hit. Three replacements per core, plus neighbor controls.
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
CORE1, CORE2 = (109, 117), (168, 176)
WT_CORE = "ACCGGATGT"
CLEAN = ["AGACGAGTG", "TCGTATCGG", "CGCATCGTA"]   # gained=0, GC-matched

seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()

EDITS = []
for lbl, (a, b) in (("core1", CORE1), ("core2", CORE2)):
    assert seq[PAD + a - 1: PAD + b] == WT_CORE
    for i, c in enumerate(CLEAN, 1):
        EDITS.append((f"{lbl}_clean{i}", a, b, c))
EDITS.append(("core1_etsmut_OLD", *CORE1, "ACCTTATGT"))   # for comparison
EDITS.append(("core2_etsmut_OLD", *CORE2, "ACCTTATGT"))

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

frames = []
for name, a, b, new in EDITS:
    ref = seq[PAD + a - 1: PAD + b]
    v = genome.Variant(chromosome="chr11", position=EL_START + a - 1,
                       reference_bases=ref, alternate_bases=new)
    print(f"scoring {name} ({ref} -> {new}) ...")
    s = model.score_variant(
        interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=v, variant_scorers=[scorer])
    df = variant_scorers.tidy_scores(s)
    df["edit"] = name
    frames.append(df)

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_clean_mutant_scores.tsv", sep="\t", index=False)

FIB = r"fibroblast|BJ|IMR-90|HFFc6|AG04450"
sub = out[out.biosample_name.astype(str).str.contains(FIB, case=False, na=False)]
piv = (sub[sub.gene_name.isin(["KLC2", "ENSG00000255320"])]
       .groupby(["edit", "gene_name"]).raw_score.mean().unstack().round(4))
pd.set_option("display.width", 200)
print("\n=== fibroblast mean score ===")
print(piv.to_string())
print("\nreference: full 217-bp deletion = +0.0339")
