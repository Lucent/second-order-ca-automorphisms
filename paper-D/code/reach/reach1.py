"""Reach of the edge-move method, part (1): no kinks next to a background with a stuck right edge.

For an elementary rule f let N_f = {(b,c): f(0,b,c) = f(1,b,c)}. If a configuration B has a stuck
right edge at r, then (Y_B(t,r+1), Y_B(t,r+2)) lies in N_f for every t.

For each rule, the program
  1. enumerates all spatially periodic configurations B of spatial period p <= PMAX with their
     F-orbits and keeps those for which some pair of adjacent columns (j, j+1) satisfies
     (Y_B(t,j), Y_B(t,j+1)) in N_f for all t (a necessary condition for a stuck right edge);
  2. for 1 <= b <= BMAX and 0 < |a| <= b builds the window graph of Fix(sigma^a F^b): vertices are
     the words of length 2b, edges the words w of length 2b+1 with F^b(w) = w_{b+a}; the graph is
     pruned to the vertices on bi-infinite paths, and its strongly connected components are computed;
  3. for every kept B all of whose windows of length 2b are vertices, checks whether a path from a
     window of B reaches a strongly connected component that contains a cycle and contains no window
     of B. Such a path is a kink of Fix(sigma^a F^b) with B on its left.
It prints, per rule, the number of kept backgrounds, the number of triples (B, a, b) tested in
step 3, and the triples with a kink (at most 5 are listed).

Usage: python reach1.py RULE [RULE ...]      (defaults PMAX = 6, BMAX = 4)
       python reach1.py classes              the 44 rules of the 11 classes of the proposition
       python reach1.py controls             rules 110, 14, 93, 25, 87, 159
"""
import sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

PMAX, BMAX = 6, 4

def bit(r, a, b, c): return (r >> ((a << 2) | (b << 1) | c)) & 1
def mk(g): return sum(g(a, b, c) << ((a << 2) | (b << 1) | c) for a in (0, 1) for b in (0, 1) for c in (0, 1))
def mirror(r): return mk(lambda a, b, c: bit(r, c, b, a))
def kappa(r): return mk(lambda a, b, c: bit(r, 1 - a, 1 - b, 1 - c))
def cls(r): return sorted({r, mirror(r), kappa(r), mirror(kappa(r))})

def cells_of(idx, n):
    return np.stack([((idx >> (2 * (n - 1 - i))) & 3) for i in range(n)], axis=1)

def step(rule, W):
    x = W >> 1; y = W & 1
    v = (rule >> (4 * y[:, :-2] + 2 * y[:, 1:-1] + y[:, 2:])) & 1
    return (y[:, 1:-1] << 1) | (x[:, 1:-1] ^ v)

def backgrounds(rule):
    """list of (p, tuple of cells) for the kept periodic configurations"""
    Nf = {(b, c) for b in (0, 1) for c in (0, 1) if bit(rule, 0, b, c) == bit(rule, 1, b, c)}
    out = []
    for p in range(1, PMAX + 1):
        idx = np.arange(4 ** p, dtype=np.int64); C = cells_of(idx, p)
        ext = np.concatenate([C[:, -1:], C, C[:, :1]], axis=1)
        FC = step(rule, ext); P = np.zeros(len(idx), dtype=np.int64)
        for i in range(p): P = (P << 2) | FC[:, i]
        y = C & 1; mask = np.zeros(len(idx), dtype=np.int64)
        for j in range(p):
            ok = np.zeros(len(idx), dtype=bool)
            for (b, c) in Nf: ok |= (y[:, j] == b) & (y[:, (j + 1) % p] == c)
            mask |= ok.astype(np.int64) << j
        Q = P.copy()
        for _ in range(2 * p + 2):          # AND of the mask over the whole F-orbit (orbit length <= 4^p)
            mask &= mask[Q]; Q = Q[Q]
        for s in np.nonzero(mask)[0]: out.append((p, tuple(int(v) for v in C[s])))
    return out

def fix_graph(rule, a, b):
    L = 2 * b + 1; idx = np.arange(4 ** L, dtype=np.int64); W = cells_of(idx, L); V = W
    for _ in range(b): V = step(rule, V)
    e = idx[V[:, 0] == W[:, b + a]]
    src = e >> 2; dst = e & (4 ** (2 * b) - 1); n = 4 ** (2 * b)
    alive = np.zeros(n, dtype=bool); alive[src] = True; alive[dst] = True
    while True:                              # prune to the vertices on bi-infinite paths
        k = alive[src] & alive[dst]; src, dst = src[k], dst[k]
        na = np.zeros(n, dtype=bool); hs = np.zeros(n, dtype=bool); hp = np.zeros(n, dtype=bool)
        hs[src] = True; hp[dst] = True; na = hs & hp
        if (na == alive).all(): break
        alive = na
    G = csr_matrix((np.ones(len(src), dtype=np.int8), (src, dst)), shape=(n, n))
    ncomp, lab = connected_components(G, directed=True, connection='strong')
    size = np.bincount(lab, minlength=ncomp)
    rec = np.zeros(ncomp, dtype=bool); rec[size > 1] = True; rec[lab[src[src == dst]]] = True
    rec[np.unique(lab[~alive])] = False
    return n, alive, src, dst, lab, rec

def can_reach(n, src, dst, target):
    """boolean vector: vertices from which some vertex of target is reachable"""
    seen = target.copy(); front = target.copy()
    while front.any():
        nxt = np.zeros(n, dtype=bool); nxt[src[front[dst]]] = True
        front = nxt & ~seen; seen |= front
    return seen

def run(rule):
    B = backgrounds(rule); tested = 0; kinks = []
    for b in range(1, BMAX + 1):
        for a in list(range(-b, 0)) + list(range(1, b + 1)):
            n, alive, src, dst, lab, rec = fix_graph(rule, a, b)
            reach = {}
            for (p, cfg) in B:
                wins = sorted({sum(cfg[(j + i) % p] << (2 * (2 * b - 1 - i)) for i in range(2 * b)) for j in range(p)})
                if not all(alive[w] for w in wins): continue
                tested += 1
                own = {int(lab[w]) for w in wins}
                for k in np.nonzero(rec)[0]:
                    k = int(k)
                    if k in own: continue
                    if k not in reach: reach[k] = can_reach(n, src, dst, lab == k)
                    if any(reach[k][w] for w in wins):
                        kinks.append((a, b, p, ''.join(map(str, cfg)))); break
    print(f"rule {rule}: {len(B)} backgrounds kept, {tested} triples (B,a,b) tested, {len(kinks)} with a kink"
          + (f"; first: {kinks[:5]}" if kinks else ""), flush=True)
    return len(kinks)

if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["classes"]:
        rules = sum((cls(r) for r in (28, 29, 41, 43, 47, 57, 61, 111, 143, 155, 157)), [])
    elif args == ["controls"]:
        rules = [110, 14, 93, 25, 87, 159]
    else:
        rules = [int(x) for x in args]
    total = sum(run(r) for r in rules)
    print("total number of triples with a kink:", total)
