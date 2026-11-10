"""GF(2) kernel dimension of the commutation system for a LINEAR second-order rule (90 or 150) at radius R,
using union-find for the 2-term equalities and dense bit-packed numpy elimination for the rest.
dim(kernel) = log2(number of block maps of radius R commuting with F). Prediction for 90R: 2^{2R}+1."""
import sys, time
import numpy as np
from centralizer import make_rule

def main(rule, R):
    t0 = time.time()
    f, F_loc, _ = make_rule(rule)
    n = 2 * R + 1; m = 2 * R + 3; N = 4 ** n; W = 4 ** m
    idx = np.arange(W, dtype=np.int64)
    cells = [((idx >> (2 * (m - 1 - i))) & 3).astype(np.int8) for i in range(m)]
    Ft = np.array([[[F_loc(a, b, c) for c in range(4)] for b in range(4)] for a in range(4)], dtype=np.int8)
    Fw = [Ft[cells[i], cells[i + 1], cells[i + 2]] for i in range(m - 2)]
    eL = sum(Fw[i].astype(np.int64) << (2 * (n - 1 - i)) for i in range(n))
    e1 = (idx >> 4) & (N - 1); e2 = (idx >> 2) & (N - 1); e3 = idx & (N - 1)
    del cells, Fw
    # linear rule: f(a,b,c) = fa*a ^ fb*b ^ fc*c (^ const). verify linearity
    fa, fb, fc, f0 = f(1,0,0)^f(0,0,0), f(0,1,0)^f(0,0,0), f(0,0,1)^f(0,0,0), f(0,0,0)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                assert f(a,b,c) == (fa&a) ^ (fb&b) ^ (fc&c) ^ f0, "rule not linear"
    assert f0 == 0
    nv = 2 * N
    parent = np.arange(nv, dtype=np.int64)
    pairs = np.unique(np.stack([2 * eL, 2 * e2 + 1], 1), axis=0)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb: parent[max(ra, rb)] = min(ra, rb)
    for _ in range(64):
        p2 = parent[parent]
        if (p2 == parent).all(): break
        parent = p2
    roots, cls = np.unique(parent, return_inverse=True); C = len(roots)
    print(f"rule {rule}R R={R}: {nv} table bits, {len(pairs)} distinct equalities, {C} classes [{time.time()-t0:.0f}s]", flush=True)
    # xor equations over classes: Y_eL ^ X_e2 ^ fa*Y_e1 ^ fb*Y_e2 ^ fc*Y_e3 = 0
    terms = [cls[2 * eL + 1], cls[2 * e2]]
    if fa: terms.append(cls[2 * e1 + 1])
    if fb: terms.append(cls[2 * e2 + 1])
    if fc: terms.append(cls[2 * e3 + 1])
    T = np.stack(terms, 1).astype(np.int64)
    # cancel repeated classes within a row (mod 2) by sorting and pairing
    T.sort(axis=1)
    words = (C + 63) // 64
    # build packed rows
    rows = np.zeros((W, words), dtype=np.uint64)
    for j in range(T.shape[1]):
        col = T[:, j]
        np.bitwise_xor.at(rows, (np.arange(W), col // 64), (np.uint64(1) << (col % 64).astype(np.uint64)))
    rows = np.unique(rows, axis=0)
    rows = rows[rows.any(axis=1)]
    print(f"  {len(rows)} distinct nonzero xor-equations over {C} classes ({words} words/row) [{time.time()-t0:.0f}s]", flush=True)
    # dense elimination, column by column
    M = rows; nrows = len(M); rank = 0; t1 = time.time()
    pivot_rows_done = np.zeros(nrows, dtype=bool)
    for col in range(C):
        w, b = col // 64, np.uint64(col % 64)
        has = ((M[:, w] >> b) & np.uint64(1)).astype(bool)
        cand = np.flatnonzero(has & ~pivot_rows_done)
        if len(cand) == 0: continue
        p = cand[0]; pivot_rows_done[p] = True; rank += 1
        has[p] = False
        tgt = np.flatnonzero(has)
        if len(tgt): M[tgt] ^= M[p]
        if col % 2048 == 0: print(f"    col {col}/{C}, rank {rank} [{time.time()-t1:.0f}s]", flush=True)
    print(f"=== rule {rule}R R={R}: classes {C}, rank {rank}, kernel dim {C - rank} (predicted for 90R: 2^(2R)+1 = {4**R+1}) [{time.time()-t0:.0f}s]", flush=True)

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
