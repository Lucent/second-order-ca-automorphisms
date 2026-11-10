#!/usr/bin/env python3
"""Check for Lemma pm:lem:theta and Proposition pm:prop:odd (holonomy behind the front).
Builds the two-copy window graph Gamma_{w,M} from the LOCAL RULE only
(F(c)_i = (y_i, x_i ^ f(y_{i-1},y_i,y_{i+1}))), takes its eventual image T,
and G_w = T with the layer y[w] (y-bit of the leftmost site) forgotten.
usage: ref_check.py RULE W M [flip] 
"""
import sys, itertools
import numpy as np

def run(rule, w, M, verbose=True):
    n = w + 1 + M                       # sites k=0..n-1 <-> site p-w+k ; front k=w
    mask = (1 << n) - 1
    tt = [(rule >> i) & 1 for i in range(8)]
    # ANF of q where f = a ^ q(b,c)
    q = lambda b, c: tt[2*b + c]
    assert all(tt[4 + 2*b + c] == 1 ^ q(b, c) for b in (0, 1) for c in (0, 1)), "not left-permutive"
    alpha = q(0, 0) ^ q(1, 0); kcoef = q(0, 0) ^ q(1, 0) ^ q(0, 1) ^ q(1, 1)

    def fbits(L, C, R):                 # bitwise f(L,C,R) on packed ints / arrays
        out = 0
        for a in (0, 1):
            for b in (0, 1):
                for c in (0, 1):
                    if tt[4*a + 2*b + c]:
                        out = out | ((L if a else ~L) & (C if b else ~C) & (R if c else ~R))
        return out

    def step1(x, y, x1, y1, y2):        # one copy; arrays of packed words
        Y = y | (y1 << n) | (y2 << (n + 1))
        X = x | (x1 << n)
        ny = ((X ^ fbits(Y << 1, Y, Y >> 1)) >> 1) & mask
        nx = (Y >> 1) & mask
        return nx, ny

    def pack(ax, ay, bx, by): return ax | (ay << n) | (bx << 2*n) | (by << 3*n)
    def unpack(s): return s & mask, (s >> n) & mask, (s >> 2*n) & mask, (s >> 3*n) & mask
    def step(S, r):
        ax, ay, bx, by = unpack(S); x1, y1, x2, y2 = r
        nax, nay = step1(ax, ay, x1, y1, y2); nbx, nby = step1(bx, by, x1, y1, y2)
        return pack(nax, nay, nbx, nby)

    # all window states: a free; b = a on k>w; at k=w same x, different y; k<w free
    lo = (1 << w) - 1
    A = np.arange(1 << (2*n), dtype=np.int64); ax = A & mask; ay = A >> n
    B = np.arange(1 << (2*w), dtype=np.int64); bxl = B & lo; byl = B >> w
    ax = np.repeat(ax, len(B)); ay = np.repeat(ay, len(B))
    bx = (ax & ~lo) | np.tile(bxl, len(A)); by = ((ay & ~lo) ^ (1 << w)) | np.tile(byl, len(A))
    S = pack(ax, ay, bx, by)
    inputs = [(a, b, 0, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]   # x(r2) is never read
    sizes = [len(S)]
    while True:
        S2 = np.unique(np.concatenate([np.unique(step(S, r)) for r in inputs]))
        sizes.append(len(S2))
        if len(S2) == len(S) and np.array_equal(S2, S): break
        S = S2
    T = S
    pred = 2**(3*w + 2) * 4**M
    # --- G_w : forget y[w] = bit k=0 of ay, by
    drop = ~((1 << n) | (1 << 3*n))
    V = np.unique(T & drop)
    assert len(V) * 4 == len(T), "y[w] not free over G_w"
    idx = {int(v): i for i, v in enumerate(V)}
    nV = len(V)
    ax, ay, bx, by = unpack(V)
    edges = {}                           # (i,(u,v)) -> (j, theta)
    for r in inputs:
        x1, y1, x2, y2 = r
        Vn = step(V, r)
        nax, nay, nbx, nby = unpack(Vn)
        # label (u,v) of Lemma r30:lem:inputs: v = x-bit... new cell at p+M+1 in copy a is (u,v)
        u = (nax >> (n - 1)) & 1; v = (nay >> (n - 1)) & 1
        phia = (nay ^ 0) & 1; phib = nby & 1          # y[w] was set to 0 in V, so new bit = increment
        # check increment independent of y[w]:
        for da, db in ((1, 0), (0, 1), (1, 1)):
            W2 = step(V | (da << n) | (db << 3*n), r)
            assert np.array_equal(W2 & drop, Vn & drop)
            assert np.array_equal(((W2 >> n) & 1) ^ da, phia) and np.array_equal(((W2 >> 3*n) & 1) ^ db, phib)
        ya = (ay >> 1) & 1; yb = (by >> 1) & 1          # layer y[w-1] = site k=1
        th = 1 | (phia << 1) | (phib << 2) | ((ya ^ alpha) << 3) | ((yb ^ alpha) << 4)
        J = Vn & drop
        for i in range(nV):
            key = (i, (int(u[i]), int(v[i])))
            val = (idx[int(J[i])], int(th[i]))
            assert edges.get(key, val) == val
            edges[key] = val
    assert len(edges) == 4 * nV
    out = [[] for _ in range(nV)]; inn = [[] for _ in range(nV)]
    for (i, lab), (j, th) in edges.items():
        out[i].append((j, th)); inn[j].append(i)
    def reach(adj):
        seen = [False]*nV; seen[0] = True; st = [0]
        while st:
            i = st.pop()
            for j in adj[i]:
                jj = j[0] if isinstance(j, tuple) else j
                if not seen[jj]: seen[jj] = True; st.append(jj)
        return sum(seen)
    sc = (reach(out) == nV and reach(inn) == nV)
    # --- exact SET of Psi-values of closed walks at vertex 0: BFS in the product G x (Z/2)^5
    seen = np.zeros((nV, 32), dtype=bool); seen[0, 0] = True; st = [(0, 0)]
    while st:
        i, g = st.pop()
        for j, th in out[i]:
            h = g ^ th
            if not seen[j, h]: seen[j, h] = True; st.append((j, h))
    vals0 = sorted(int(g) for g in np.nonzero(seen[0])[0])
    res = dict(rule=rule, w=w, M=M, alpha=alpha, k=kcoef, sizes=sizes, T=len(T), pred=pred,
               nV=nV, sc=sc, nvals=len(vals0), vals=vals0)
    hol = sorted({(g >> 1) & 3 for g in vals0}); res['hol'] = hol
    if verbose:
        print(f"rule {rule} k={kcoef} alpha={alpha} w={w} M={M}: |W|..|T| {sizes} pred {pred} ok={len(T)==pred}; "
              f"|G_w|={nV} strongly_conn={sc}; #Psi values of closed walks at a vertex = {len(vals0)}; "
              f"Hol(l_a,l_b) values={hol}", flush=True)
    return res, (n, mask, V, idx, edges, alpha, unpack, pack)

def flipcheck(rule, w, M):
    res, (n, mask, V, idx, edges, alpha, unpack, pack) = run(rule, w, M, verbose=False)
    nV = len(V)
    ax, ay, bx, by = unpack(V)
    for c in 'ab':
        # top flip: complement y_c[w-1] (y-bit of site k=1) and x_c[w-1] (x-bit of site p-w, k=0)
        if c == 'a': fm = (1 << 0) | (1 << (n + 1))
        else: fm = (1 << (2*n)) | (1 << (3*n + 1))
        img = [idx.get(int(v) ^ fm) for v in V]
        if any(i is None for i in img):
            print(f"  rule {rule} w={w} M={M} flip_{c}: does NOT preserve the vertex set "
                  f"({sum(i is None for i in img)} of {nV} images outside)"); continue
        okedge = True; dphi_ok = True; other_ok = True; dN_ok = True
        ylow = ((ay if c == 'a' else by) >> 2) & 1 if w >= 2 else None   # y_c[w-2] = site k=2
        for (i, lab), (j, th) in edges.items():
            j2, th2 = edges[(img[i], lab)]
            if j2 != img[j]: okedge = False
            d = th ^ th2
            bitc, bito, nc, no = (1, 2, 3, 4) if c == 'a' else (2, 1, 4, 3)
            if ((d >> bitc) & 1) != (int(ylow[i]) ^ alpha): dphi_ok = False
            if ((d >> bito) & 1) or ((d >> no) & 1) or (d & 1): other_ok = False
            if ((d >> nc) & 1) != 1: dN_ok = False
        print(f"  rule {rule} w={w} M={M} flip_{c}: vertex set preserved; edges commute with same label={okedge}; "
              f"phi_c changes by y_c[d-2]+alpha={dphi_ok}; N_c entry changes by 1={dN_ok}; other entries fixed={other_ok}")
    # swap
    sw = [idx.get(int(pack(bx[i], by[i], ax[i], ay[i]))) for i in range(nV)]
    ok = all(s is not None for s in sw)
    if ok:
        for (i, lab), (j, th) in edges.items():
            lab2 = (lab[0], lab[1] ^ (1 if M == 0 else 0))
            j2, th2 = edges[(sw[i], lab2)]
            exp = (th & 1) | (((th >> 2) & 1) << 1) | (((th >> 1) & 1) << 2) | (((th >> 4) & 1) << 3) | (((th >> 3) & 1) << 4)
            if j2 != sw[j] or th2 != exp: ok = False
    print(f"  rule {rule} w={w} M={M} swap: automorphism with label (u,v+[M=0]) and psi entries exchanged = {ok}")

if __name__ == '__main__':
    rule, w, M = map(int, sys.argv[1:4])
    if len(sys.argv) > 4: flipcheck(rule, w, M)
    else: run(rule, w, M)
