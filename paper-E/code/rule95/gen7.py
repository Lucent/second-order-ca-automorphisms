"""Tables of h_* = sigma^{-1} hatH_{phi_*,id} (phi_*: exchange of 0001001 and 0000101), built from the definition by
   simulation:
   TY[w]: y'_i as a function of the 15 bits  b_m = y_{i+m} (m odd), x_{i+m} (m even), |m|<=7, w = sum b_m << (m+7);
   TX[w]: x'_i as a function of              b_m = x_{i+m} (m odd), y_{i+m} (m even).
   Every sample is checked for consistency (the output is a function of these 15 bits); every entry is hit."""
import numpy as np, sys, os, random
from core import *

def debruijn(n):
    a = [0] * (2 * n); seq = []
    def db(t, p):
        if t > n:
            if n % p == 0: seq.extend(a[1:p + 1])
        else:
            a[t] = a[t - p]; db(t + 1, p)
            for j in range(a[t - p] + 1, 2):
                a[t] = j; db(t + 1, t)
    db(1, 1); return seq

def build(phi, psi, nsig, R, par_y, chunk=900, seed=0, passes=2, log=True):
    """generic: tables for hatH_{phi,psi} sigma^nsig as a radius-R sublattice map.
       par_y: y'_i reads y_{i+m} for m = par_y mod 2 and x_{i+m} otherwise; x'_i the complementary bits."""
    nb = 2 * R + 1
    T = [np.full(1 << nb, -1, dtype=np.int8), np.full(1 << nb, -1, dtype=np.int8)]   # T[0]=TX, T[1]=TY
    dB = debruijn(nb); N = len(dB); rnd = random.Random(seed); incons = 0; nsamp = 0
    offA = rnd.randrange(N); offB = rnd.randrange(N)
    for ps in range(passes):
        for start in range(0, N, chunk - 2 * nb):
            L = chunk
            # class sequences: A_k = y_k (k even), x_k (k odd);  B_k = x_k (k even), y_k (k odd)   (k = absolute cell)
            base = ps                                    # shift the parity between passes
            cells = []
            for i in range(L):
                k = base + i
                A = dB[(offA + start + i) % N]; B = dB[(offB + start + i) % N]
                cells.append((B, A) if k % 2 == 0 else (A, B))
            c = from_cells(base, cells)
            lo2, X2, Y2 = hatH(*c, phi, psi, nsig=nsig)
            img = to_dict(lo2, X2, Y2); src = to_dict(*c)
            get = lambda d, k, b: d.get(k, (1, 1))[b]
            for k in range(base - 2 * nb, base + L + 2 * nb):
                for out in (0, 1):                       # 0: x', 1: y'
                    w = 0
                    for m in range(-R, R + 1):
                        same = (m - par_y) % 2 == 0      # y' reads y at m = par_y (mod 2)
                        layer = (1 if same else 0) if out == 1 else (0 if same else 1)
                        w |= get(src, k + m, layer) << (m + R)
                    v = get(img, k, out); nsamp += 1
                    if T[out][w] == -1: T[out][w] = v
                    elif T[out][w] != v: incons += 1
        if log: print("pass", ps, "unfilled:", int((T[0] < 0).sum()), int((T[1] < 0).sum()), "inconsistent samples:", incons, "of", nsamp, flush=True)
    assert (T[0] >= 0).all() and (T[1] >= 0).all(), 'table not covered'
    return T[0], T[1], incons

if __name__ == "__main__":
    phi1 = gm_swap('0001001', '0000101')
    TX, TY, inc = build(phi1, ID, nsig=-1, R=7, par_y=1)
    assert inc == 0 and (TX >= 0).all() and (TY >= 0).all()
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))     # output directory
    np.save(os.path.join(out, "tx7.npy"), TX.astype(np.uint8)); np.save(os.path.join(out, "ty7.npy"), TY.astype(np.uint8))
    TX.astype(np.uint8).tofile(os.path.join(out, "tx7.bin")); TY.astype(np.uint8).tofile(os.path.join(out, "ty7.bin"))
    print("saved; ones:", int(TX.sum()), int(TY.sum()))
