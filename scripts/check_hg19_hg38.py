#!/usr/bin/env python3
"""Locate the published hg19 interval's sequence within hg38. Stdlib only."""
import json, subprocess

def fetch(genome, start0, end):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};"
           f"chrom=chr11;start={start0};end={end}")
    return json.loads(subprocess.run(["curl","-sSL",url],
           capture_output=True, text=True).stdout)["dna"].upper()

hg19 = fetch("hg19", 66024556, 66024773)          # published interval
win  = fetch("hg38", 66255000, 66259000)          # generous hg38 window
rc   = hg19[::-1].translate(str.maketrans("ACGT","TGCA"))

print(f"hg19 {len(hg19)} bp, ends {hg19[0]}/{hg19[-1]}")
for label, q in (("forward", hg19), ("revcomp", rc)):
    i = win.find(q)
    if i >= 0:
        s = 66255001 + i
        print(f"  {label}: FULL 217/217 match at hg38 chr11:{s}-{s+216}")
    else:
        j = win.find(q[:30])
        print(f"  {label}: no full match; first 30 bp "
              + (f"at hg38 chr11:{66255001+j}" if j >= 0 else "absent"))
print("\ncurrently assumed element: chr11:66,256,955-66,257,171")

