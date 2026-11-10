"""Rule 110R: finite checks around the proposition "no walls in 110R". All three are exhaustive.

(1) For every value of the four cells c[r-3..r], the first t in {1,2,3} at which F^t(c)_r depends on
    the cells to the right of r. Expected: t = 1, 2, 3 for 192, 48 and 16 of the 256 values.
(2) For every window c[j-T-1..j+T+1] and every change of the cell j (three changes), whether the
    two configurations differ only in the column j in all rows -1 <= t <= T. Expected: 147456 pairs
    for T = 4 and none for T = 5.
(3) For P = 1..6: the number of words c[r-P..r] such that F^P(c)_r = c_r for all values of the cells
    c[r+1..r+P]. (F^P(c)_r depends only on c[r-P..r+P], so a word of any length whose right end cell
    is fixed by F^P in every context gives such a word.) Expected: 0 for every P.

Usage: python nowalls110.py
"""
import numpy as np

RULE = 110
T8 = np.array([(RULE >> k) & 1 for k in range(8)], dtype=np.uint8)

def words(n, lo=0, hi=None):
    hi = 4 ** n if hi is None else hi
    idx = np.arange(lo, hi, dtype=np.int64)
    return np.stack([((idx >> (2 * (n - 1 - i))) & 3).astype(np.uint8) for i in range(n)], axis=1)

def step(W):
    x = W >> 1; y = W & 1
    v = T8[4 * y[:, :-2] + 2 * y[:, 1:-1] + y[:, 2:]]
    return (y[:, 1:-1] << 1) | (x[:, 1:-1] ^ v)

def power(W, t):
    for _ in range(t): W = step(W)
    return W

# ---- (1)
first = {1: 0, 2: 0, 3: 0, None: 0}
for left in range(256):
    lw = [(left >> (2 * (3 - i))) & 3 for i in range(4)]
    ft = None
    for t in (1, 2, 3):
        # window c[r-t..r+t]: the last t+1 cells of the left word, then all right contexts of length t
        R = words(t)
        W = np.concatenate([np.tile(np.array(lw[3 - t:], dtype=np.uint8), (len(R), 1)), R], axis=1)
        if len(np.unique(power(W, t)[:, 0])) > 1: ft = t; break
    first[ft] += 1
print("(1) first t at which the right context matters:", first)

# ---- (2)
def confined(W, D, T):
    """rows of (W, D) that agree outside the centre column in all rows t <= T (row -1 and 0 are the cells)"""
    ok = np.ones(len(W), dtype=bool)
    for t in range(1, T + 1):
        W = step(W); D = step(D); c = W.shape[1] // 2
        diff = W != D; diff[:, c] = False
        # rows -1..t at the columns kept: a cell (x, y) at time t holds rows t-1 and t
        ok &= ~diff.any(axis=1)
    return ok

def count(T, base=None):
    n = 2 * T + 3; j = T + 1; total = 0; keep = []
    if base is None:
        chunks = [(lo, min(lo + 4 ** 9, 4 ** n)) for lo in range(0, 4 ** n, 4 ** 9)]
        for lo, hi in chunks:
            W = words(n, lo, hi)
            for d in (1, 2, 3):
                D = W.copy(); D[:, j] ^= d
                ok = confined(W, D, T); total += int(ok.sum()); keep.append((W[ok], d))
    else:
        E = words(2)
        for W0, d in base:
            if len(W0) == 0: continue
            W = np.concatenate([np.repeat(E[:, :1], len(W0), axis=0), np.tile(W0, (16, 1)),
                                np.repeat(E[:, 1:], len(W0), axis=0)], axis=1)
            D = W.copy(); D[:, j] ^= d
            ok = confined(W, D, T); total += int(ok.sum()); keep.append((W[ok], d))
    return total, keep

n4, keep4 = count(4)
print("(2) pairs whose difference stays in its column through t = 4:", n4)
n5, _ = count(5, keep4)      # a pair confined through t = 5 restricts to one confined through t = 4
print("(2) pairs whose difference stays in its column through t = 5:", n5)

# ---- (3)
for P in range(1, 7):
    nl = 4 ** (P + 1); R = words(P); found = 0
    for lo in range(0, nl, 64):
        L = words(P + 1, lo, min(lo + 64, nl))
        W = np.concatenate([np.repeat(L, len(R), axis=0), np.tile(R, (len(L), 1))], axis=1)
        out = power(W, P)[:, 0].reshape(len(L), len(R))
        found += int((out == L[:, -1:]).all(axis=1).sum())
    print(f"(3) P = {P}: words c[r-P..r] whose right end cell is fixed by F^P in every right context: {found}")
