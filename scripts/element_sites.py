#!/usr/bin/env python3
"""Map FIMO hits in the corrected element by position cluster and TF family."""

import re, pandas as pd

FA = "raw/element_TRUE_flank500_hg38.fa"
PAD, N = 500, 217
seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
el = seq[PAD:PAD + N]

core = "ACCGGATGT"
rc = core[::-1].translate(str.maketrans("ACGT", "TGCA"))
print("ETS core occurrences (1-based in element):")
for label, q in (("+", core), ("-", rc)):
    for m in re.finditer(f"(?={q})", el):
        p = m.start()
        print(f"  {label} {p+1}-{p+len(q)}   context: {el[max(0,p-12):p]}[{q}]{el[p+len(q):p+len(q)+12]}")

d = pd.read_csv("processed/fimo_TRUE/fimo.tsv", sep="\t", comment="#").dropna(subset=["motif_alt_id"])
d = d[(d.start >= PAD + 1) & (d.stop <= PAD + N)].copy()
d["el_start"] = d.start - PAD
d["el_stop"] = d.stop - PAD

FAM = [("ETS", r"^(ETS|ETV|ELK|ERF|FEV|FLI|GABPA|EHF|ELF|SPI)"),
       ("KRAB-ZNF", r"^ZNF|^Zfp|^ZKSCAN|^ZSCAN"),
       ("NFKB", r"^NFKB|^REL"),
       ("ZBTB", r"^ZBTB|^Zbtb"),
       ("INSM", r"^INSM"),
       ("CTCF", r"^CTCF")]
def fam(m):
    for name, pat in FAM:
        if re.search(pat, m, re.I):
            return name
    return "other"
d["family"] = d.motif_alt_id.map(fam)

print("\nhits by family (collapsed by family+position):")
u = d.drop_duplicates(["family", "el_start", "el_stop"])
print(u.family.value_counts().to_string())

print("\nfamily occupancy by position:")
pd.set_option("display.width", 200)
print(u.sort_values("el_start")[["family", "motif_alt_id", "el_start", "el_stop",
                                 "strand", "p-value"]].to_string(index=False))

