import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

EL_S, EL_E = 66256955, 66257171
MOT_S, MOT_E = 66257035, 66257065
TSS = 66257270

peaks = [
    ("cortical-interneuron (KCl-)", 66256969, 66257149, "neural"),
    ("cortical-interneuron (KCl+)", 66256964, 66257090, "neural"),
    ("DND41 (T-ALL)",               66256926, 66257247, "other"),
    ("hESC_D5",                     66256971, 66257489, "other"),
    ("OCI-Ly1 (B-cell lymphoma)",   66257014, 66257503, "other"),
]
peaks.sort(key=lambda p: p[2] - p[1])

fig, ax = plt.subplots(figsize=(9, 5))
XMIN, XMAX = 66256850, 66257560

ax.axvspan(MOT_S, MOT_E, color="#FFD166", alpha=.55, zorder=0)
ax.axvspan(EL_S, EL_E, facecolor="none", edgecolor="#118AB2",
           lw=1.5, ls="--", zorder=1)

y = len(peaks) + 1.2
ax.barh(y, EL_E - EL_S, left=EL_S, height=.45,
        color="#118AB2", ec="black", lw=.7)
ax.text(EL_S - 12, y, "216-bp element", ha="right", va="center", fontsize=9)
ax.barh(y, MOT_E - MOT_S, left=MOT_S, height=.45,
        color="#EF476F", ec="black", lw=.7)

for i, (lab, s, e, kind) in enumerate(peaks):
    mid = (s + e) // 2
    col = "#06D6A0" if kind == "neural" else "#B0B0B0"
    ax.barh(i, e - s, left=s, height=.45, color=col, ec="black", lw=.7)
    ax.plot([mid], [i], marker="v", color="black", ms=7, zorder=3)
    ax.text(XMIN - 12, i, f"{lab}\n{e-s} bp", ha="right", va="center", fontsize=8)
    ax.text(e + 15, i, f"mid {mid:,}", va="center", fontsize=7, color="#444")

ax.axvline(TSS, color="black", lw=1, ls=":")
ax.text(TSS + 8, -1.0, "KLC2 TSS", fontsize=8, rotation=90, va="bottom")

ax.set_xlim(XMIN, XMAX); ax.set_ylim(-1.2, y + .9)
ax.set_yticks([]); ax.set_xlabel("chr11 (hg38)")
for sp in ["left", "right", "top"]:
    ax.spines[sp].set_visible(False)
ax.set_title("CTCF motif and ReMap ChIP-seq peaks at the SPOAN element",
             fontsize=11)

ax.legend(handles=[
    mpatches.Patch(color="#118AB2", label="216-bp deleted element"),
    mpatches.Patch(color="#EF476F", label="CTCF motif (FIMO p=1.6e-5)"),
    mpatches.Patch(color="#06D6A0", label="ReMap peak — neural"),
    mpatches.Patch(color="#B0B0B0", label="ReMap peak — non-neural"),
], fontsize=8, loc="lower right", framealpha=.95)

ax.ticklabel_format(axis="x", style="plain", useOffset=False)
ax.set_xticks([66256900, 66257000, 66257100, 66257200, 66257300, 66257400, 66257500])
ax.set_xticklabels([f"{t:,}" for t in ax.get_xticks()], fontsize=8)

plt.tight_layout()
plt.savefig("figures/ctcf_motif_peaks.png", dpi=300)
plt.savefig("figures/ctcf_motif_peaks.svg")
print("saved")

