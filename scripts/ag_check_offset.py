#!/usr/bin/env python3
"""Verify the element offset in the flank FASTA using the known CTCF motif."""

FLANK_FA = "raw/klc2_del_flank500_hg38.fa"
EL_START = 66256955
MOTIF = "CTCCAAATACCTTAACTGACACAAGGGAGCA"   # chr11:66,257,035-66,257,065 hg38
MOTIF_START = 66257035

seq = "".join(l.strip() for l in open(FLANK_FA) if not l.startswith(">")).upper()
print(f"flank file length: {len(seq)}  (expect 1217 for 500+217+500)")

i = seq.find(MOTIF)
if i < 0:
    from re import findall
    rc = seq[::-1].translate(str.maketrans("ACGT", "TGCA"))
    print("motif NOT found on this strand; revcomp hit:", rc.find(MOTIF))
else:
    pad = i - (MOTIF_START - EL_START)
    print(f"motif found at index {i}  =>  element starts at index {pad}")
    el = seq[pad:pad + 217]
    print(f"element first/last base: {el[0]} / {el[-1]}")
    print(f"5' 20: {el[:20]}")
    print(f"3' 20: {el[-20:]}")

