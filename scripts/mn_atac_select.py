#!/usr/bin/env python3
"""Reduce passing files to one pooled file per experiment, then map to donors."""

import pandas as pd

df = pd.read_csv("processed/mn_atac_passing_files.tsv", sep="\t")
df["bio_reps"] = df["bio_reps"].fillna("").astype(str)
df["n_reps"] = df["bio_reps"].apply(lambda s: len([x for x in s.split(",") if x]))

# pooled file = the one covering the most biological replicates
sel = (df.sort_values(["experiment", "n_reps"], ascending=[True, False])
         .groupby("experiment", as_index=False)
         .first())

sel.to_csv("processed/mn_atac_selected.tsv", sep="\t", index=False)

pd.set_option("display.width", 220)
print(sel[["experiment", "file", "bio_reps", "donors", "description"]].to_string(index=False))
print()
print("experiments:", sel.experiment.nunique())
print("donors:     ", sel.donors.nunique())
print()
print("experiments per donor:")
print(sel.groupby("donors").experiment.count().to_string())
