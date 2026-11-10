#!/usr/bin/env python3
"""(a) per-component order of the Psi group for an affine rule; (b) the period-3
history below 2 for rule 135 as a closed walk of G_2 (and its Psi), for several M and w;
(c) lifting: for 135, closed walks of G_{d-1} with ell=(0,0) lift to closed walks of same length.
usage: ref_base.py percomp RULE W M | base | lift W M"""
import sys, random
from ref_check import run

def graph(rule, w, M):
    res, (n, mask, V, idx, edges, alpha, unpack, pack) = run(rule, w, M, verbose=False)
    return res, n, V, idx, edges, alpha, pack

if sys.argv[1] == 'percomp':
    rule, w, M = map(int, sys.argv[2:5])
    res, n, V, idx, edges, alpha, pack = graph(rule, w, M)
    nV = len(V); out = [[] for _ in range(nV)]
    for (i, lab), (j, th) in edges.items(): out[i].append((j, th))
    import numpy as np
    orders = {}
    done = [False]*nV
    for s in range(nV):
        if done[s]: continue
        seen = np.zeros((nV, 32), bool); seen[s, 0] = True; st = [(s, 0)]
        while st:
            i, g = st.pop()
            for j, th in out[i]:
                if not seen[j, g ^ th]: seen[j, g ^ th] = True; st.append((j, g ^ th))
        vals = frozenset(int(x) for x in np.nonzero(seen[s])[0])
        orders[len(vals)] = orders.get(len(vals), 0) + 1; done[s] = True
        if len(orders) and sum(orders.values()) >= 64: break
    print(f"rule {rule} w={w} M={M}: number of Psi-values of closed walks at a vertex -> #vertices sampled: {orders}")

elif sys.argv[1] == 'base':
    q = lambda b, c: 1 ^ (b & c)            # rule 135
    Za = {-2: '101', -1: '110', 0: '011', 1: '010'}; Zb = {-2: '101', -1: '110', 0: '100', 1: '111'}
    tail = {0: '011', 2: '110', 1: '101'}   # layers i = 0,-1,-2 mod 3  (i mod 3 = 0,2,1)
    def Z(c, t, i):
        d = Za if c == 'a' else Zb
        s = d[i] if i in d else tail[i % 3]
        return int(s[t % 3])
    # (Z1), (Z2) by hand formula
    for c in 'ab':
        for t in range(3):
            for i in range(-12, 2):
                assert Z(c, t+1, i) == Z(c, t, i) ^ Z(c, t-1, i-2) ^ q(Z(c, t, i-1), Z(c, t, i-2)), (c, t, i)
    assert all(Z('a', t, i) == Z('b', t, i) for t in range(3) for i in range(-12, 0)) and all(Z('a', t, 0) != Z('b', t, 0) for t in range(3))
    print("history below 2 (rule 135, tail 011): (Z1) for i<=1 and (Z2) hold")
    for w, M in ((2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (3, 1), (4, 0)):
        res, n, V, idx, edges, alpha, pack = graph(135, w, M)
        def eta(t):
            w_ = {}
            for c in 'ab':
                x = y = 0
                for k in range(n):
                    d = w - k                       # site p-d: cell (x[d-1], y[d])
                    if d <= 1: y |= Z(c, t, d) << k
                    if d - 1 <= 1: x |= Z(c, t - 1, d - 1) << k
                w_[c] = (x, y)
            return w_
        # project to the layers of G_2: for w>2 the layers y[2..w-1], x[2..w-1] and y[w] are not part of G_2;
        # for w=2 the vertex of G_w is eta with y[2] dropped (set to 0)
        if w == 2:
            ids = []
            for t in range(4):
                e = eta(t); code = int(pack(e['a'][0], e['a'][1], e['b'][0], e['b'][1]))
                ids.append(idx.get(code))
            ok = all(i is not None for i in ids) and ids[3] == ids[0]
            th = 0
            for t in range(3):
                nxt = [(j, thv) for (i, lab), (j, thv) in edges.items() if i == ids[t] and j == ids[t+1]]
                ok = ok and len(nxt) == 1
                if nxt: th ^= nxt[0][1]
            print(f"  w={w} M={M}: closed 3-walk in G_2: {ok}; Psi_2 = (T,l_a,l_b,N_a,N_b) = {[(th>>b)&1 for b in range(5)]}")
        else:
            # w>2: check that the G_2-layers of the walk occur as the projection of a closed... (S_nu does not depend on w):
            lowmask_bits = [k for k in range(n) if w - k <= 1]
            def proj(code_x, code_y):
                x = sum(((code_x >> k) & 1) << k for k in range(n) if w - k - 1 <= 1)
                y = sum(((code_y >> k) & 1) << k for k in lowmask_bits)
                return x, y
            import numpy as np
            m = (1 << n) - 1
            PV = set()
            for v in V:
                v = int(v); PV.add(proj(v & m, (v >> n) & m) + proj((v >> 2*n) & m, (v >> 3*n) & m))
            ok = all((eta(t)['a'] + eta(t)['b']) in PV for t in range(3))
            print(f"  w={w} M={M}: the three vertices lie in the projection of S to the layers of G_2: {ok}; |proj|={len(PV)} (expect {2**(3*2+2-2)*4**M})")

elif sys.argv[1] == 'lift':
    w, M = int(sys.argv[2]), int(sys.argv[3])
    res, n, V, idx, edges, alpha, pack = graph(135, w, M)
    # G_w -> G_{w-1}: forget y[w-1] (y-bit k=1) and x[w-1] (x-bit k=0) in both copies
    fm = (1 << 0) | (1 << (n + 1)) | (1 << 2*n) | (1 << (3*n + 1))
    nV = len(V); out = {}
    for (i, lab), (j, th) in edges.items(): out[(i, lab)] = (j, th)
    base = [int(v) & ~fm for v in V]
    # check: projected edge target depends only on the projection (so G_{w-1} is the quotient)
    q = {}
    for (i, lab), (j, th) in edges.items():
        key = (base[i], lab); assert q.setdefault(key, base[j]) == base[j]
    # phi at step y[w-1] on G_{w-1}: increment of y[w-1] (k=1 y-bits) along an edge
    def phi_prev(i, j):
        a = int(V[i]) ^ int(V[j]); return ((a >> (n + 1)) & 1, (a >> (3*n + 1)) & 1)
    rnd = random.Random(1); labs = [(u, v) for u in (0, 1) for v in (0, 1)]
    tested = good = 0
    for trial in range(20000):
        i0 = rnd.randrange(nV); L = rnd.randrange(1, 9); i = i0; ell = [0, 0]; seq = []
        for _ in range(L):
            lab = rnd.choice(labs); j, th = out[(i, lab)]; p = phi_prev(i, j); ell[0] ^= p[0]; ell[1] ^= p[1]; seq.append(lab); i = j
        if base[i] != base[i0]: continue          # projection not closed
        if ell != [0, 0]: 
            # then the lift must NOT close in the y[w-1] layer
            continue
        tested += 1
        # lift from zeta: run the label sequence once (arbitrary x[w-1]) then again
        j = i
        for lab in seq: j = out[(j, lab)][0]
        k = j
        for lab in seq: k = out[(k, lab)][0]
        good += (k == j)
    print(f"lift check 135 w={w} M={M}: closed walks of G_{w-1} with ell=(0,0) sampled: {tested}; lifted closed walk of the same length through zeta: {good}")
