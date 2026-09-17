#!/usr/bin/env python3
"""Which genes move, within lineage? Averaging across opposite-signed
tissue classes cancels the effect and is uninterpretable."""

import pandas as pd
d = pd.read_csv("processed/ag_klc2_scores.tsv", sep="\t")
d = d[d.variant == "del217"]

CLASSES = {
    "fibroblast": r"fibroblast|BJ|IMR-90|HFFc6|AG04450",
    "nerve_cord": r"tibial nerve|sciatic nerve|spinal cord|motor neuron",
    "lymphoid":   r"T cell|B cell|natural killer|lymphoblast|GM12878|PBMC|"
                  r"peripheral blood mononuclear|venous blood|thymus|spleen",
}

pd.set_option("display.width", 200)
for name, pat in CLASSES.items():
    sub = d[d.biosample_name.astype(str).str.contains(pat, case=False, na=False)]
    g = (sub.groupby("gene_name").raw_score
            .agg(mean="mean", n="size").reset_index())
    g["abs"] = g["mean"].abs()
    print(f"\n=== {name}  ({sub.biosample_name.nunique()} biosamples, "
          f"{len(sub)} rows) — top 12 ===")
    print(g.sort_values("abs", ascending=False).head(12).round(4).to_string(index=False))
    rank = g.sort_values("abs", ascending=False).reset_index(drop=True)
    for gene in ["KLC2", "ENSG00000255320"]:
        r = rank.index[rank.gene_name == gene]
        print(f"  {gene}: rank {r[0]+1 if len(r) else 'absent'} of {len(rank)}")

