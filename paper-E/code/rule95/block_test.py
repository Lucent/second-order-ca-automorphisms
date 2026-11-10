"""F_95 has no blocking word: explicit witnesses for every word w of length <= LMAX and every position q.
   c = w in the background (1,1); c' = c plus one bare defect of Phi^+ (y_{-a}=0, x_{-a-1}=0) far to the left.
   Both contain w at the same place; the program finds the first n >= 0 with F^n(c)_q != F^n(c')_q.
   Also: the same with a pseudo-random pair of contexts.   usage: block_test.py [LMAX]"""
import numpy as np, sys
from core import FT
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
def step(X, Y):
    Yl = np.concatenate((np.ones((Y.shape[0], 1), np.uint8), Y[:, :-1]), axis=1); Yr = np.concatenate((Y[:, 1:], np.ones((Y.shape[0], 1), np.uint8)), axis=1)
    return Y, X ^ FT[Yl, Y, Yr]
rng = np.random.default_rng(5)
for L in range(1, LMAX + 1):
    nw = 4 ** L; D = 4 * L + 12; T = 12 * L + 4 * D + 40; W = L + 2 * T + 2 * D + 20; o = W // 2
    idx = np.arange(nw)
    X = np.ones((nw, W), np.uint8); Y = np.ones((nw, W), np.uint8)
    for i in range(L):
        s = (idx >> (2 * (L - 1 - i))) & 3; X[:, o + i] = s >> 1; Y[:, o + i] = s & 1
    worst = 0; fails = 0
    for dpar in (0, 1):
        X2, Y2 = X.copy(), Y.copy(); pos = o - D - dpar           # cell -a
        Y2[:, pos] = 0; X2[:, pos - 1] = 0
        A, B = (X.copy(), Y.copy()), (X2, Y2)
        first = np.full((nw, L), -1)
        for n in range(T):
            diff = (A[0][:, o:o + L] != B[0][:, o:o + L]) | (A[1][:, o:o + L] != B[1][:, o:o + L])
            first = np.where((first < 0) & diff, n, first)
            A = step(*A); B = step(*B)
        if dpar == 0: res = first
        else: res = np.where(res < 0, first, np.where(first < 0, res, np.minimum(res, first)))
    fails = int((res < 0).sum()); worst = int(res.max())
    # random contexts
    Xa = rng.integers(0, 2, (nw, W), dtype=np.uint8); Ya = rng.integers(0, 2, (nw, W), dtype=np.uint8)
    Xb = rng.integers(0, 2, (nw, W), dtype=np.uint8); Yb = rng.integers(0, 2, (nw, W), dtype=np.uint8)
    for M in (Xa, Xb): M[:, o:o + L] = X[:, o:o + L]; M[:, :T + 2] = 1; M[:, -T - 2:] = 1
    for M in (Ya, Yb): M[:, o:o + L] = Y[:, o:o + L]; M[:, :T + 2] = 1; M[:, -T - 2:] = 1
    A, B = (Xa, Ya), (Xb, Yb); first = np.full((nw, L), -1)
    for n in range(60):
        diff = (A[0][:, o:o + L] != B[0][:, o:o + L]) | (A[1][:, o:o + L] != B[1][:, o:o + L])
        first = np.where((first < 0) & diff, n, first); A = step(*A); B = step(*B)
    print("length %d: %d words x %d positions; defect witness missing for %d (word,position) pairs, largest first time of difference %d; random pair of contexts: missing %d, largest first time %d" % (L, nw, L, fails, worst, int((first < 0).sum()), int(first.max())), flush=True)
