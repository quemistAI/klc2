import pyBigWig, statistics

EL_S, EL_E = 66024556, 66024773        # hg19 BED
L_S, L_E   = EL_S-200, EL_S            # 200 bp left flank
R_S, R_E   = EL_E, EL_E+200            # 200 bp right flank

for eid in ["E073","E055","E029"]:
    try:
        bw = pyBigWig.open(f"raw/signal/{eid}-H3K27ac.fc.signal.bigwig")
    except Exception as e:
        print(f"{eid}: no file ({e})"); continue
    def mean(s,e):
        v = [x for x in bw.values("chr11", s, e) if x is not None]
        return statistics.mean(v) if v else float("nan")
    el, l, r = mean(EL_S,EL_E), mean(L_S,L_E), mean(R_S,R_E)
    flank = (l+r)/2
    print(f"{eid}: element {el:.2f} | left {l:.2f} | right {r:.2f} "
          f"| flank mean {flank:.2f} | ratio {el/flank:.2f}")
    bw.close()
