"""Exact verification: tau o tau = id and tau o D != D o tau, by composing the radius-3 tables on ALL windows of
length 13 (4^13 = 67M) via numpy; also tau o F = F o tau re-checked and the 24 changed entries listed."""
import numpy as np, json, time
from centralizer import make_rule, index_word
from classify_r4 import check_commutes_np_R
f, F_loc, Finv_loc = make_rule(54); R = 3; n = 7; N = 4 ** n
PI = {0: 1, 1: 2, 2: 0, 3: 3}
def tau_local(w):
    c = w[3]
    if c == 3 and w[1] == 3 and w[2] == 3 and w[4] != 3 and w[5] == 3 and w[6] == 3: return PI[w[4]]
    if c != 3 and w[1] == 3 and w[2] == 3 and w[4] != 3 and w[5] == 3 and w[6] == 3 and c == PI[w[4]]: return 3
    return c
def D_local(w): return PI[w[3]] if (w[1] == 3 and w[2] == 3 and w[3] != 3 and w[4] == 3 and w[5] == 3) else w[3]
T = np.array([tau_local(index_word(e, n)) for e in range(N)], dtype=np.int64)
D = np.array([D_local(index_word(e, n)) for e in range(N)], dtype=np.int64)
I = np.array([index_word(e, n)[3] for e in range(N)], dtype=np.int64)
def compose13(g, h):
    """g o h on all 13-windows, returned as an array over 4^13 (centre = position 6)."""
    out = np.empty(4 ** 13, dtype=np.int8); t0 = time.time()
    for start in range(0, 4 ** 13, 1 << 22):
        idx = np.arange(start, min(4 ** 13, start + (1 << 22)), dtype=np.int64)
        hs = []
        for p in range(7):                     # h at positions 3..9 uses cells p..p+6 of the 13-window
            hs.append(h[(idx >> (2 * (12 - (p + 6)))) & (N - 1)])
        gi = np.zeros(len(idx), dtype=np.int64)
        for v in hs: gi = (gi << 2) | v
        out[start:start + len(idx)] = g[gi]
    return out
centre13 = lambda: ((np.arange(4 ** 13, dtype=np.int64) >> 12) & 3).astype(np.int8)
TT = compose13(T, T); print("tau o tau == id on all 13-windows:", bool((TT == centre13()).all()), flush=True)
TD = compose13(T, D); DT = compose13(D, T); print("tau o D == D o tau:", bool((TD == DT).all()), "; windows where they differ:", int((TD != DT).sum()), flush=True)
print("tau commutes with F:", check_commutes_np_R(list(T), F_loc, R), "; D commutes with F:", check_commutes_np_R(list(D), F_loc, R))
changed = [(''.join(map(str, index_word(e, n))), int(T[e])) for e in range(N) if T[e] != I[e]]
print(len(changed), "changed entries:", changed)
