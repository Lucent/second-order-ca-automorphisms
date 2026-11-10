"""For every particle-conserving block map of G of window 7 and 8 that the graph test declares non-injective,
   exhibit two distinct periodic points of G (period <= PMAX) with the same image.   usage: inj_witness.py [PMAX]"""
import sys, numpy as np
from enum_golden import gwords, injective, enumerate_maps
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
def cyc(n): return [m for m in range(1 << n) if not (m & (m >> 1)) and not ((m & 1) and (m >> (n - 1)) & 1)]
for w in (7, 8):
    sols, _ = enumerate_maps(w); E = gwords(w)
    lut = np.zeros((len(sols), 1 << w), dtype=np.uint8)
    for k, f in enumerate(sols):
        for e in E: lut[k, e] = f[e]
    inj = np.array([injective(w, f) for f in sols])
    found = np.zeros(len(sols), dtype=bool); per = np.zeros(len(sols), dtype=int)
    for n in range(1, PMAX + 1):
        pts = np.array(cyc(n), dtype=np.int64)
        bits = [(pts >> (n - 1 - i)) & 1 for i in range(n)]
        win = []
        for i in range(n):
            e = np.zeros(len(pts), dtype=np.int64)
            for k in range(w): e = (e << 1) | bits[(i + k) % n]
            win.append(e)
        win = np.stack(win, axis=1)                        # points x n
        for k in np.nonzero(~found)[0]:
            img = lut[k][win]                              # points x n
            code = img.astype(np.int64) @ (1 << np.arange(n, dtype=np.int64))
            if len(np.unique(code)) < len(pts): found[k] = True; per[k] = n
    print("window %d: %d maps, injective by the graph test %d; non-injective %d, with an explicit collision of periodic points %d (largest period needed %d); injective maps with a collision: %d" % (
        w, len(sols), inj.sum(), (~inj).sum(), (found & ~inj).sum(), per[~inj].max(), (found & inj).sum()), flush=True)
