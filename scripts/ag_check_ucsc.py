#!/usr/bin/env python3
"""Verify the flank FASTA's absolute genomic registration against UCSC."""

import requests

FLANK_FA = "raw/klc2_del_flank500_hg38.fa"
EL_START, EL_END = 66256955, 66257171     # 1-based inclusive

seq = "".join(l.strip() for l in open(FLANK_FA) if not l.startswith(">")).upper()
local = seq[500:717]

def ucsc(start1, end1):   # 1-based inclusive in, UCSC 0-based half-open out
    r = requests.get("https://api.genome.ucsc.edu/getSequence",
                     params={"genome": "hg38", "chrom": "chr11",
                             "start": start1 - 1, "end": end1}, timeout=30)
    r.raise_for_status()
    return r.json()["dna"].upper()

for shift in (-1, 0, 1):
    ref = ucsc(EL_START + shift, EL_END + shift)
    tag = "MATCH" if ref == local else "     "
    print(f"shift {shift:+d}: {ref[:20]} ... {ref[-20:]}  first/last {ref[0]}/{ref[-1]}  {tag}")

print(f"\nlocal      : {local[:20]} ... {local[-20:]}  first/last {local[0]}/{local[-1]}")

# context around both boundaries
print("\n5' context 66,256,950-66,256,960:", ucsc(66256950, 66256960))
print("3' context 66,257,166-66,257,176:", ucsc(66257166, 66257176))
