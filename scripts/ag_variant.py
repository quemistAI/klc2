#!/usr/bin/env python3
"""AlphaGenome calibration test: does the SPOAN deletion raise predicted KLC2?
Runs both the 217-bp and 216-bp interpretations of the deletion.
RUN IN: conda env `alphagenome`.  Reads raw/, writes processed/."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ag_env import check_env

API_KEY = check_env()

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

EL_START, EL_END = 66256955, 66257171          # hg38, 1-based inclusive, 217 bp
FLANK_FA = "raw/klc2_del_flank500_hg38.fa"     # element +/- 500 bp
FLANK_PAD = 500
ONTOLOGY = [
    # paste the exact IDs from scripts/ag_ontology.py output
    "UBERON:0001323",   # tibial nerve            <- VERIFY
    "UBERON:0000178",   # whole blood             <- VERIFY
]

# --- load reference sequence -------------------------------------------------
seq = "".join(l.strip() for l in open(FLANK_FA) if not l.startswith(">")).upper()
el = seq[FLANK_PAD:FLANK_PAD + (EL_END - EL_START + 1)]
assert len(el) == 217, f"element is {len(el)} bp, expected 217"
print(f"element 5' 20bp: {el[:20]}")
print(f"element 3' 20bp: {el[-20:]}")
print(f"first base {el[0]}, last base {el[-1]} "
      f"(both G => 216/217 ambiguity applies)")

anchor = seq[FLANK_PAD - 1]

VARIANTS = {
    "del217": genome.Variant(
        chromosome="chr11", position=EL_START - 1,
        reference_bases=anchor + el, alternate_bases=anchor),
    "del216_drop3p": genome.Variant(
        chromosome="chr11", position=EL_START - 1,
        reference_bases=anchor + el[:-1], alternate_bases=anchor),
}

# --- predict -----------------------------------------------------------------
model = dna_client.create(API_KEY)
scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

frames = []
for label, v in VARIANTS.items():
    print(f"\n=== {label}: {len(v.reference_bases)-1} bp deleted ===")
    interval = v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    scores = model.score_variant(interval=interval, variant=v,
                                 variant_scorers=[scorer])
    try:
        df = variant_scorers.tidy_scores(scores)
    except AttributeError:
        # fallback: unpack AnnData manually
        import numpy as np
        parts = []
        for s in scores:
            m = pd.DataFrame(s.X, index=s.obs.index, columns=s.var.index)
            m = m.stack().rename("raw_score").reset_index()
            m.columns = ["gene_row", "track_row", "raw_score"]
            for c in s.obs.columns:
                m[c] = s.obs.loc[m.gene_row, c].values
            for c in s.var.columns:
                m[c] = s.var.loc[m.track_row, c].values
            parts.append(m)
        df = pd.concat(parts, ignore_index=True)
    df["variant"] = label
    frames.append(df)

out = pd.concat(frames, ignore_index=True)
out.to_csv("processed/ag_klc2_scores.tsv", sep="\t", index=False)

print("\ncolumns:", list(out.columns))
print("shape:", out.shape)

gene_col = next((c for c in out.columns if "gene_name" in c.lower()), None)
if gene_col:
    klc2 = out[out[gene_col].astype(str).str.fullmatch("KLC2", case=False, na=False)]
    keep = [c for c in ["variant", gene_col, "biosample_name", "ontology_curie",
                        "output_type", "raw_score", "quantile_score"]
            if c in klc2.columns]
    pd.set_option("display.width", 200)
    pd.set_option("display.max_rows", 200)
    print("\n--- KLC2 ---")
    print(klc2[keep].sort_values("raw_score", ascending=False).to_string(index=False))
else:
    print("\nno gene_name column; here are the non-numeric columns:")
    print([c for c in out.columns if out[c].dtype == object])
