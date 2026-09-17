#!/usr/bin/env python3
"""Fetch an hg38 interval (+/- 500 bp flanks) from UCSC as FASTA.
Usage: python3 scripts/fetch_hg38.py chr11 66254800 66255016 raw/ctrl1.fa
Coordinates are 1-based inclusive, as displayed in the UCSC browser.
Stdlib only - runs in any conda env."""

import sys, json, subprocess

PAD = 500
chrom, start, end, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
core_len = end - start + 1

url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};"
       f"start={start - 1 - PAD};end={end + PAD}")
raw = subprocess.run(["curl", "-sS", url], capture_output=True, text=True).stdout
try:
    dna = json.loads(raw)["dna"].upper()
except Exception:
    sys.exit(f"UCSC returned non-JSON:\n{raw[:300]}")

expect = core_len + 2 * PAD
assert len(dna) == expect, f"got {len(dna)} bp, expected {expect}"

core = dna[PAD:PAD + core_len]
gc = 100 * (core.count("G") + core.count("C")) / len(core)

with open(out, "w") as f:
    f.write(f">{chrom}:{start}-{end}_hg38_flank{PAD} core={core_len}bp\n")
    for i in range(0, len(dna), 60):
        f.write(dna[i:i + 60] + "\n")

print(f"{out}")
print(f"  core {core_len} bp, GC {gc:.1f}%  (element = 54.4%)")
print(f"  core 5' {core[:15]} ... 3' {core[-15:]}")
raw = subprocess.run(["curl", "-sSL", url], capture_output=True, text=True).stdout
