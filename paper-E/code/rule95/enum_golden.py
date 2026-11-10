"""All block maps phi: G -> G of the golden mean shift with window w that conserve the number of 1s of finite points.
   (phi x)_i = f(x_i ... x_{i+w-1}).  Exact criterion: phi conserves iff there is lam: {G-words of length w-1} -> Z,
   lam(0..0) = 0, with  f(e) - e_0 = lam(e_0..e_{w-2}) - lam(e_1..e_{w-1})  for every G-word e of length w
   (the weight f(e)-e_0 sums to 0 on every closed walk through 0^{w-1} iff it does on every cycle of the strongly
   connected de Bruijn graph of G iff it is a coboundary).  So the maps correspond bijectively to the integer
   functions lam with lam(0)=0, e_0 + lam(src) - lam(tgt) in {0,1} on every edge, and f(e)f(e')=0 for consecutive e,e'.
   Injectivity: exact, by pruning the graph of pairs of paths with equal images."""
import sys, itertools
from collections import deque

def gwords(n): return [m for m in range(1 << n) if not (m & (m >> 1))]

def enumerate_maps(w):
    V = gwords(w - 1); E = gwords(w); mask = (1 << (w - 1)) - 1
    src = {e: e >> 1 for e in E}; tgt = {e: e & mask for e in E}; e0 = {e: (e >> (w - 1)) & 1 for e in E}
    out = {v: [] for v in V}; inn = {v: [] for v in V}
    for e in E: out[src[e]].append(e); inn[tgt[e]].append(e)
    # BFS order from 0 using edges in both directions
    order = [0]; seen = {0}; dq = deque([0])
    while dq:
        v = dq.popleft()
        for e in out[v] + inn[v]:
            for u in (src[e], tgt[e]):
                if u not in seen: seen.add(u); order.append(u); dq.append(u)
    assert len(order) == len(V)
    pos = {v: i for i, v in enumerate(order)}
    # constraints attached to the last-assigned vertex involved
    edge_at = [[] for _ in V]
    for e in E: edge_at[max(pos[src[e]], pos[tgt[e]])].append(e)
    pair_at = [[] for _ in V]                         # consecutive edges e -> e' (a G-word of length w+1)
    for e in E:
        for e2 in out[tgt[e]]:
            pair_at[max(pos[src[e]], pos[tgt[e]], pos[tgt[e2]])].append((e, e2))
    lam = {}; sols = []; nodes = [0]
    def fval(e): return e0[e] + lam[src[e]] - lam[tgt[e]]
    def rec(k):
        nodes[0] += 1
        if k == len(order):
            sols.append({e: fval(e) for e in E}); return
        v = order[k]
        # candidate values from the edges to vertices already assigned
        cands = None
        for e in out[v]:
            if tgt[e] in lam and tgt[e] != v:
                c = {lam[tgt[e]] - e0[e], lam[tgt[e]] - e0[e] + 1}; cands = c if cands is None else cands & c
        for e in inn[v]:
            if src[e] in lam and src[e] != v:
                c = {lam[src[e]] + e0[e], lam[src[e]] + e0[e] - 1}; cands = c if cands is None else cands & c
        for val in sorted(cands):
            lam[v] = val
            if all(fval(e) in (0, 1) for e in edge_at[k]) and all(not (fval(a) and fval(b)) for a, b in pair_at[k]):
                rec(k + 1)
            del lam[v]
    lam[0] = 0
    if all(fval(e) in (0, 1) for e in edge_at[0]) and all(not (fval(a) and fval(b)) for a, b in pair_at[0]): rec(1)
    return sols, nodes[0]

def injective(w, f):
    """exact: phi is injective on G iff the pair graph has no off-diagonal vertex lying on a bi-infinite path."""
    V = gwords(w - 1); E = gwords(w); mask = (1 << (w - 1)) - 1
    out = {v: [] for v in V}
    for e in E: out[e >> 1].append(e)
    nodes = [(a, b) for a in V for b in V]
    succ = {p: [] for p in nodes}; pred = {p: [] for p in nodes}
    for (a, b) in nodes:
        for e in out[a]:
            for e2 in out[b]:
                if f[e] == f[e2]:
                    q = (e & mask, e2 & mask); succ[(a, b)].append(q); pred[q].append((a, b))
    alive = set(nodes)
    changed = True
    while changed:
        changed = False
        for p in list(alive):
            if not any(q in alive for q in succ[p]) or not any(q in alive for q in pred[p]):
                alive.discard(p); changed = True
    return not any(a != b for (a, b) in alive)

def describe(w, f):
    """identify shifts"""
    E = gwords(w)
    for k in range(w):
        if all(f[e] == (e >> (w - 1 - k)) & 1 for e in E): return "sigma^%d" % k
    return None

def table_str(w, f): return "".join(str(f[e]) for e in gwords(w))

if __name__ == "__main__":
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for w in range(2, W + 1):
        sols, nodes = enumerate_maps(w)
        inj = [f for f in sols if injective(w, f)]
        names = [describe(w, f) for f in inj]
        print("w=%d: G-words %d, conserving maps %d (search nodes %d), injective %d: %s" % (w, len(gwords(w)), len(sols), nodes, len(inj), names), flush=True)
        for f in inj:
            if describe(w, f) is None:
                diff = [format(e, '0%db' % w) + "->" + str(f[e]) for e in gwords(w) if f[e] != (e >> (w - 1 - 4)) & 1] if w > 4 else []
                print("   non-shift injective map; entries differing from sigma^4:", diff)
