import pyBigWig
import pandas as pd
import matplotlib.pyplot as plt

EIDS = {
    "E073": "Neural - DLPFC",
    "E055": "Fibroblast",
    "E029": "Myeloid - monocyte",
    "E030": "Myeloid - neutrophil",
    "E124": "Myeloid - CD14+ monocyte",
    "E032": "Lymphoid - B cell",
    "E034": "Lymphoid - T cell",
    "E047": "Lymphoid - CD8+ naive",
}
MARKS = ["H3K4me3", "H3K4me1", "H3K27me3"]

CHROM = "chr11"
EL_START, EL_END = 66024556, 66024773      # hg19, BED coords
PAD = 1500
WIN_START, WIN_END = EL_START - PAD, EL_END + PAD

rows = []
for eid, label in EIDS.items():
    for mark in MARKS:
        path = f"raw/signal/{eid}-{mark}.fc.signal.bigwig"
        try:
            bw = pyBigWig.open(path)
            vals = bw.values(CHROM, WIN_START, WIN_END)
            bw.close()
        except Exception as e:
            print(f"FAILED {eid} {mark}: {e}")
            continue
        for i, v in enumerate(vals):
            rows.append({"eid": eid, "label": label, "mark": mark,
                         "pos": WIN_START + i,
                         "val": 0.0 if v is None else float(v)})
        print(f"ok {eid} {mark}")

df = pd.DataFrame(rows)
df.to_csv("processed/histone_signal.csv", index=False)

# mean signal over the element vs the rest of the window
el = df[(df.pos >= EL_START) & (df.pos < EL_END)]
summary = (el.groupby(["eid", "label", "mark"])["val"]
             .mean().round(3).reset_index())
print("\nMean fold-change over the 217-bp element:")
print(summary.to_string(index=False))
summary.to_csv("processed/histone_element_means.csv", index=False)

fig, axes = plt.subplots(len(EIDS), 1, figsize=(10, 18), sharex=True)
for ax, (eid, label) in zip(axes, EIDS.items()):
    sub = df[df.eid == eid]
    for mark in MARKS:
        m = sub[sub.mark == mark]
        ax.plot(m.pos, m.val, label=mark, lw=1)
    ax.axvspan(EL_START, EL_END, color="cyan", alpha=0.25)
    ax.set_ylabel("fold change")
    ax.set_title(f"{eid} — {label}", fontsize=9, loc="left")
    ax.legend(fontsize=7, ncol=3)
axes[-1].set_xlabel(f"{CHROM} (hg19)")
plt.tight_layout()
plt.savefig("figures/histone_signal.png", dpi=300)
plt.savefig("figures/histone_signal.svg")
print("\ndone")
