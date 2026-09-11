import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

df = pd.read_csv("processed/chromhmm/klc2_states.tsv", sep="\t")

meta = {
    "E073": ("Neural", "E073 DLPFC"),
    "E081": ("Neural", "E081 Fetal brain"),
    "E053": ("Neural", "E053 Neurosphere ctx"),
    "E054": ("Neural", "E054 Neurosphere GE"),
    "E055": ("Fibroblast", "E055 Fibroblast"),
    "E056": ("Fibroblast", "E056 Fibroblast"),
    "E032": ("Lymphoid", "E032 B cell"),
    "E034": ("Lymphoid", "E034 T cell"),
    "E046": ("Lymphoid", "E046 NK cell"),
    "E047": ("Lymphoid", "E047 CD8+ naive"),
    "E035": ("Stem", "E035 HSC"),
    "E116": ("Lymphoid", "E116 GM12878*"),
    "E062": ("Myeloid-mix", "E062 PBMC"),
    "E029": ("Myeloid", "E029 Monocyte"),
    "E030": ("Myeloid", "E030 Neutrophil"),
    "E124": ("Myeloid", "E124 CD14+ monocyte"),
}

# collapse multiple overlapping segments per EID, keeping the
# most repressive state observed
priority = ["13_ReprPC", "14_ReprPCWk", "12_EnhBiv", "11_BivFlnk",
            "10_TssBiv", "9_Het", "15_Quies", "7_Enh", "6_EnhG",
            "2_TssAFlnk", "1_TssA"]
rank = {s: i for i, s in enumerate(priority)}
df["rank"] = df["state"].map(lambda s: rank.get(s, 99))
best = df.sort_values("rank").groupby("EID").first().reset_index()

best["group"] = best["EID"].map(lambda e: meta[e][0])
best["label"] = best["EID"].map(lambda e: meta[e][1])

order = ["Neural", "Fibroblast", "Stem", "Lymphoid",
         "Myeloid-mix", "Myeloid"]
best["gorder"] = best["group"].map({g: i for i, g in enumerate(order)})
best = best.sort_values(["gorder", "EID"])

colors = {
    "1_TssA": "#FF0000", "2_TssAFlnk": "#FF4500",
    "11_BivFlnk": "#E9967A", "12_EnhBiv": "#BDB76B",
    "13_ReprPC": "#808080", "14_ReprPCWk": "#C0C0C0",
    "15_Quies": "#FFFFFF",
}

fig, ax = plt.subplots(figsize=(5.5, 7))
for i, (_, r) in enumerate(best.iterrows()):
    ax.barh(i, 1, color=colors.get(r["state"], "#DDDDDD"),
            edgecolor="black", linewidth=0.6)
    ax.text(1.03, i, r["state"], va="center", fontsize=8)

ax.set_yticks(range(len(best)))
ax.set_yticklabels(best["label"], fontsize=8)
ax.invert_yaxis()
ax.set_xticks([]); ax.set_xlim(0, 1.6)
ax.set_title("ChromHMM state at the SPOAN element\n"
             "(chr11:66,024,557-66,024,773, hg19)", fontsize=10)

handles = [mpatches.Patch(color=c, label=s, ec="black")
           for s, c in colors.items() if s in set(best["state"])]
ax.legend(handles=handles, fontsize=7, loc="lower right",
          title="State", title_fontsize=7)

plt.tight_layout()
plt.savefig("figures/chromhmm_heatmap.png", dpi=300)
plt.savefig("figures/chromhmm_heatmap.svg")
print(best[["label", "group", "state"]].to_string(index=False))
