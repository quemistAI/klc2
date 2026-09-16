import pyBigWig, statistics, glob, os

EL   = (66256954, 66257171)
LEFT = (66256754, 66256954)
PROM = (66257300, 66257600)
BKG  = (66250000, 66251000)

EXCLUDE = {"ENCFF638TCZ"}   # different normalisation: background 0.01
                            # vs ~0.56; promoter 81.6 vs ~2.9

def mean(bw, s, e):
    v = [x for x in bw.values("chr11", s, e) if x is not None]
    return statistics.mean(v) if v else float("nan")

rows = []
for path in sorted(glob.glob("raw/mn_atac/*.bigWig")):
    fid = os.path.basename(path).split(".")[0]
    bw = pyBigWig.open(path)
    el, lf, pr, bg = (mean(bw, *r) for r in (EL, LEFT, PROM, BKG))
    bw.close()
    tag = "  [EXCLUDED]" if fid in EXCLUDE else ""
    print(f"{fid}  element {el:7.2f} | left {lf:7.2f} | promoter {pr:7.2f} "
          f"| background {bg:7.2f} | el/bkg {el/bg:7.2f} | el/prom {el/pr:5.2f}{tag}")
    if fid not in EXCLUDE:
        rows.append((el, pr, bg))

if rows:
    eb = [r[0]/r[2] for r in rows]
    ep = [r[0]/r[1] for r in rows]
    print(f"\nincluded files: {len(rows)}")
    print(f"  element/background  {statistics.mean(eb):.2f}"
          f"  (range {min(eb):.2f}-{max(eb):.2f})")
    print(f"  element/promoter    {statistics.mean(ep):.2f}"
          f"  (range {min(ep):.2f}-{max(ep):.2f})")
