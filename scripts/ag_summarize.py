#!/usr/bin/env python3
"""Aggregate AlphaGenome KLC2 scores per biosample; check locus specificity."""

import pandas as pd
d = pd.read_csv("processed/ag_klc2_scores.tsv", sep="\t")

k = d[d.gene_name.astype(str).str.fullmatch("KLC2", case=False, na=False)]
agg = (k.groupby(["variant", "biosample_name"])
         .agg(n_tracks=("raw_score", "size"),
              mean=("raw_score", "mean"),
              lo=("raw_score", "min"),
              hi=("raw_score", "max"),
              q_mean=("quantile_score", "mean"))
         .round(4).reset_index())
agg.to_csv("processed/ag_klc2_by_biosample.tsv", sep="\t", index=False)

FOCUS = ("fibroblast|BJ|IMR-90|tibial nerve|sciatic|spinal cord|motor neuron|"
         "cerebellum|blood|T cell|B cell|lymphoblast|GM12878|monocyte|PBMC|"
         "peripheral blood|spleen|thymus|natural killer")
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 300)
print("=== disease-relevant + blood ===")
print(agg[agg.biosample_name.str.contains(FOCUS, case=False, na=False)]
      .sort_values(["variant", "mean"], ascending=[True, False]).to_string(index=False))

# is the effect specific to KLC2, or does everything in the window move?
print("\n=== top 20 genes by |mean score|, del217 ===")
g = (d[d.variant == "del217"].groupby("gene_name").raw_score
       .agg(["mean", "size"]).reset_index())
g["abs"] = g["mean"].abs()
print(g.sort_values("abs", ascending=False).head(20).round(4).to_string(index=False))
