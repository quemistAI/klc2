#!/usr/bin/env python3
"""Base-resolution ATAC profile across the SPOAN element and KLC2 promoter.
Remote reads - no bigWig downloads. Decides whether the element is a discrete
feature or the shoulder of the promoter peak."""

import pyBigWig, pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CHR, START, END, NBINS = "chr11", 66255600, 66258600, 300
EL = (66257086, 66257302)

urls = pd.read_csv("processed/mn_atac_urls.tsv", sep="\t")
x = np.linspace(START, END, NBINS)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
traces = {}

for _, r in urls.iterrows():
    print(f"  {r.file} ...", flush=True)
    bw = pyBigWig.open(r.url)
    v = np.nan_to_num(np.array(bw.stats(CHR, START, END, nBins=NBINS, type="mean"),
                               dtype=float))
    bw.close()
    traces[r.file] = v
    ax1.plot(x, v, lw=0.9, alpha=0.6, label=f"{r.donor}")

# per-file max-normalised, so donors with different depth are comparable
for f, v in traces.items():
    ax2.plot(x, v / (v.max() or 1), lw=0.9, alpha=0.6)

for ax in (ax1, ax2):
    ax.axvspan(*EL, color="grey", alpha=0.3, zorder=0)
ax1.set_ylabel("ATAC fold change over control")
ax2.set_ylabel("normalised to per-file max")
ax2.set_xlabel("chr11 position (hg38)")
ax1.legend(fontsize=6, ncol=3, title="donor")
ax1.set_title("iPSC-derived motor neuron ATAC across the SPOAN element (10 exp, 6 donors)")

fig.tight_layout()
fig.savefig("figures/mn_atac_profile.png", dpi=300)
np.savez_compressed("processed/mn_atac_profile.npz", x=x, **traces)
print("wrote figures/mn_atac_profile.png")
