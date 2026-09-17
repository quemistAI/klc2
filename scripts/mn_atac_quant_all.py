#!/usr/bin/env python3
"""Quantify ATAC signal over the SPOAN element by remote range reads.
No bigWig downloads. Falls back to bigWigSummary if pyBigWig lacks curl."""

import pandas as pd, numpy as np, subprocess, sys

CHR = "chr11"
WINDOWS = {                          # hg38, corrected 2026-09-17
    "background": (66250000, 66251000),
    "left":       (66256586, 66257085),
    "element":    (66257086, 66257302),
    "promoter":   (66257303, 66257802),
}
BKG_FLOOR = 0.05   # prespecified 2026-09-15

try:
    import pyBigWig
    USE_PYBW = bool(pyBigWig.remote)
except ImportError:
    USE_PYBW = False
print("method:", "pyBigWig (remote)" if USE_PYBW else "bigWigSummary", file=sys.stderr)


def mean_signal(url, chrom, start, end):
    if USE_PYBW:
        bw = pyBigWig.open(url)
        try:
            v = bw.stats(chrom, start, end, type="mean")[0]
        finally:
            bw.close()
        return float(v or 0.0)
    out = subprocess.run(
        ["bigWigSummary", "-type=mean", url, chrom, str(start), str(end), "1"],
        capture_output=True, text=True)
    t = out.stdout.strip()
    return float(t) if t and t not in ("n/a", "") else 0.0


urls = pd.read_csv("processed/mn_atac_urls.tsv", sep="\t")

rows = []
for _, r in urls.iterrows():
    print(f"  {r.file} ...", file=sys.stderr, flush=True)
    try:
        v = {k: mean_signal(r.url, CHR, s, e) for k, (s, e) in WINDOWS.items()}
    except Exception as e:
        print(f"    FAILED: {e}", file=sys.stderr)
        continue
    rows.append({
        "experiment": r.experiment, "file": r.file, "donor": r.donor, **v,
        "el_bkg":  v["element"] / v["background"] if v["background"] else np.nan,
        "el_left": v["element"] / v["left"]       if v["left"]       else np.nan,
        "el_prom": v["element"] / v["promoter"]   if v["promoter"]   else np.nan,
        "passes_floor": v["background"] >= BKG_FLOOR,
    })

df = pd.DataFrame(rows).round(3)
df.to_csv("processed/mn_atac_quant_all_TRUE.tsv", sep="\t", index=False)

pd.set_option("display.width", 220)
print(df.to_string(index=False))

ok = df[df.passes_floor]
print(f"\nexcluded by background floor (<{BKG_FLOOR}): {(~df.passes_floor).sum()}")

per_donor = ok.groupby("donor")[["el_bkg", "el_left", "el_prom"]].mean()
print("\nper donor:")
print(per_donor.round(3).to_string())

print(f"\nn files {len(ok)} | n experiments {ok.experiment.nunique()} | n donors {len(per_donor)}")
for col in ["el_bkg", "el_left", "el_prom"]:
    s = per_donor[col]
    print(f"{col:8s} median {s.median():.2f}  range {s.min():.2f}-{s.max():.2f}")

