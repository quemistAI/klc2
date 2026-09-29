#!/usr/bin/env python3
"""Reconstruct the SPOAN deletion junction at corrected coordinates."""

FA, PAD, N = "raw/element_TRUE_flank500_hg38.fa", 500, 217
seq = "".join(l.strip() for l in open(FA) if not l.startswith(">")).upper()
left, right = seq[:PAD], seq[PAD + N:]
junction = left[-100:] + right[:100]

el = seq[PAD:PAD + N]
assert len(junction) == 200
assert el[:30] not in junction and el[-30:] not in junction, "element leaked in"

with open("raw/junction_TRUE_200bp.fa", "w") as f:
    f.write(">KLC2_SPOAN_junction_TRUE_100bp_each_side\n")
    for i in range(0, 200, 60):
        f.write(junction[i:i+60] + "\n")
print("junction:", junction)
