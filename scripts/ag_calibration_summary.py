#!/usr/bin/env python3
"""Tissue pattern and gene specificity from the corrected-coordinate run."""

import pandas as pd
d = pd.read_csv("processed/ag_calibration_TRUE.tsv", sep="\t")
pd.set_option("display.width", 200)

CLASSES = {
    "fibroblast": r"fibroblast|BJ|IMR-90|HFFc6|AG04450",
    "nerve_cord": r"tibial nerve|sciatic nerve|spinal cord",
    "motor_neuron": r"motor neuron",
    "cerebellum": r"cerebellum",
    "lymphoid": r"T cell|B cell|natural killer|lymphoblast|GM12878|"
                r"peripheral blood mononuclear|venous blood|thymus|spleen",
    "myeloid": r"monocyte|neutrophil|macrophage",
}

# --- KLC2 by tissue class, del217 ---
k = d[(d.gene_name.astype(str).str.fullmatch("KLC2", na=False)) & (d.variant == "del217")]
print("=== KLC2 by tissue class (del217) ===")
rows = []
for name, pat in CLASSES.items():
    sub = k[k.biosample_name.astype(str).str.contains(pat, case=False, na=False)]
    if len(sub):
        per_bio = sub.groupby("biosample_name").raw_score.mean()
        rows.append({"class": name, "n_biosamples": len(per_bio),
                     "mean": per_bio.mean().round(4),
                     "lo": per_bio.min().round(4), "hi": per_bio.max().round(4),
                     "q_mean": sub.quantile_score.mean().round(4)})
print(pd.DataFrame(rows).to_string(index=False))

# --- gene specificity within lineage ---
print("\n=== gene rank within lineage (del217) ===")
for name in ["fibroblast", "nerve_cord", "lymphoid"]:
    sub = d[(d.variant == "del217") &
            (d.biosample_name.astype(str).str.contains(CLASSES[name], case=False, na=False))]
    g = sub.groupby("gene_name").raw_score.mean().reset_index()
    g["abs"] = g.raw_score.abs()
    g = g.sort_values("abs", ascending=False).reset_index(drop=True)
    r = g.index[g.gene_name == "KLC2"]
    print(f"\n{name}: {len(g)} genes | KLC2 rank "
          f"{r[0]+1 if len(r) else 'absent'}")
    print(g.head(5).round(4).to_string(index=False))

# --- 216 vs 217 agreement ---
print("\n=== del216 vs del217, KLC2 fibroblast ===")
fib = d[(d.gene_name.astype(str).str.fullmatch("KLC2", na=False)) &
        (d.biosample_name.astype(str).str.contains(CLASSES["fibroblast"], case=False, na=False))]
print(fib.groupby("variant").raw_score.mean().round(4).to_string())
