import pandas as pd

df = pd.read_csv("processed/fimo/fimo.tsv", sep="\t", comment="#")
df = df.dropna(subset=["motif_alt_id"])
print(f"total hits: {len(df)}   unique TFs: {df.motif_alt_id.nunique()}")

priority = [
    "EZH2","SUZ12","MTF2","JARID2","RNF2",
    "REST","ZBTB7A","ZBTB33","E2F1","E2F4","E2F6",
    "SP1","SP2","SP3","KLF","EGR1","MAZ","ZNF143",
    "YY1","CTCF","NRF1","GABPA",
]
mask = df.motif_alt_id.str.upper().str.contains("|".join(priority), na=False)
hits = df[mask].sort_values("p-value")

cols = ["motif_alt_id","start","stop","strand","score","p-value","matched_sequence"]
print("\nPrioritized candidates:")
print(hits[cols].to_string(index=False))
hits.to_csv("processed/fimo/fimo_prioritized.csv", index=False)

print("\nTop 20 by p-value regardless of prior:")
print(df.sort_values("p-value").head(20)[cols].to_string(index=False))
