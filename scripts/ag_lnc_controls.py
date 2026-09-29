#!/usr/bin/env python3
"""Control: matched 217-bp deletions elsewhere in the ENSG00000255320 intron.
Tests whether the predicted lncRNA effect is specific to the SPOAN element or
an artifact of deleting within the gene body.
RUN IN: conda env `alphagenome`."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env
API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

PAD = 500

DELS = {
    "element": ("raw/element_TRUE_flank500_hg38.fa", 66257086, 66257302),
    "ctrl1":   ("raw/lnc_ctrl1_flank500_hg38.fa",    66245849, 66246065),
    "ctrl2":   ("raw/lnc_ctrl2_flank500_hg38.fa",    66251399, 66251615),
    "ctrl3":   ("raw/lnc_ctrl3_flank500_hg38.fa",    66253949, 66254165),
}

model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

frames = []
for label, (fa, start, end) in DELS.items():
    if not os.path.exists(fa):
        print(f"MISSING {fa} - run the fetch_hg38.py command for {label}")
        continue
    seq = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
    n = end - start + 1
    core, anchor = seq[PAD:PAD + n], seq[PAD - 1]
    assert len(core) == n, f"{label}: core is {len(core)} bp, expected {n}"
    gc = 100 * (core.count("G") + core.count("C")) / n
    v = genome.Variant(chromosome="chr11", position=start - 1,
                       reference_bases=anchor + core, alternate_bases=anchor)
    print(f"scoring {label}: {n} bp at chr11:{start}-{end}, GC {gc:.1f}% ...")
    s = model.score_variant(
        interval=v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB),
        variant=v, variant_scorers=[scorer])
    df = variant_scorers.tidy_scores(s)
    df["deletion"] = label
    df["gc"] = round(gc, 1)
    frames.append(df)

if not frames:
    sys.exit("nothing scored - check the FASTA files exist in raw/")

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_lnc_control_scores.tsv", sep="\t", index=False)

FIB = r"fibroblast|BJ|IMR-90|HFFc6|AG04450"
sub = out[out.biosample_name.astype(str).str.contains(FIB, case=False, na=False)]
piv = (sub[sub.gene_name.isin(["KLC2", "ENSG00000255320"])]
       .groupby(["deletion", "gene_name"]).raw_score.mean().unstack().round(4))
pd.set_option("display.width", 200)
print("\n=== fibroblast mean score ===")
print(piv.to_string())
print("\nreference: element gave lncRNA -0.1303, KLC2 +0.0339")
