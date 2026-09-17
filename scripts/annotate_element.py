#!/usr/bin/env python3
"""Recompute basic locus facts at the corrected element coordinates."""

EL_START, EL_END = 66257086, 66257302
CPG_START, CPG_END = 66257440, 66258782
OLD_START, OLD_END = 66256955, 66257171
CTCF = (66257035, 66257065)

print(f"element          chr11:{EL_START:,}-{EL_END:,}  ({EL_END-EL_START+1} bp)")
print(f"previous (wrong) chr11:{OLD_START:,}-{OLD_END:,}")
ov = max(0, min(EL_END, OLD_END) - max(EL_START, OLD_START) + 1)
print(f"  offset {EL_START-OLD_START} bp, overlap {ov} bp "
      f"({100*ov/(EL_END-EL_START+1):.0f}%)\n")

print(f"CpG island       chr11:{CPG_START:,}-{CPG_END:,}")
print(f"  island starts {CPG_START-EL_END-1} bp downstream of element end")
print(f"  (was {CPG_START-OLD_END-1} bp under the old coordinates)\n")

print(f"CTCF motif       chr11:{CTCF[0]:,}-{CTCF[1]:,}")
print(f"  ends {EL_START-CTCF[1]-1} bp BEFORE element start - outside the deletion\n")

print("TO FILL IN FROM UCSC at chr11:66,256,600-66,257,800 (hg38):")
for q in ["distance to MANE KLC2 TSS",
          "distance to nearest annotated TSS",
          "ENCODE cCRE overlap (or gap boundaries)",
          "position within ENSG00000255320 intron",
          "Multiz: mouse/dog/elephant alignment present?",
          "RepeatMasker overlap"]:
    print(f"  [ ] {q}")
