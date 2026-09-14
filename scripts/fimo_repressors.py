import pandas as pd

df = pd.read_csv("processed/fimo/fimo.tsv", sep="\t", comment="#")
df = df.dropna(subset=["motif_alt_id"])
cols = ["motif_alt_id","start","stop","strand","score","p-value","q-value","matched_sequence"]

# KRAB zinc fingers: ZNF/ZKSCAN/ZBTB/ZFP families are predominantly repressive
krab = df[df.motif_alt_id.str.upper().str.match(r"(ZNF|ZKSCAN|ZBTB|ZFP|ZSCAN)", na=False)]
print("=== KRAB / C2H2 zinc fingers ===")
print(krab.sort_values("p-value")[cols].to_string(index=False) if len(krab) else "  none")

# named repressors / corepressor-recruiting factors
named = ["REST","NRSF","CTCF","YY1","SNAI","ZEB","TWIST","HIC1","PRDM",
         "BCL6","MECP2","NCOR","SIN3","EZH","SUZ12","MTF2","JARID",
         "E2F4","E2F6","THRB","NR2F","RXR","MXI1","MAD","MNT"]
mask = df.motif_alt_id.str.upper().str.contains("|".join(named), na=False)
print("\n=== Named repressors / corepressor recruiters ===")
print(df[mask].sort_values("p-value")[cols].to_string(index=False))

print(f"\nTotal hits: {len(df)} | KRAB-family: {len(krab)} | named repressors: {mask.sum()}")
