"""Reach of the edge-move method, part (1), in the form of the proposition: uniqueness of the
continuation of a background with a stuck right edge inside Fix(sigma^a F^b).

For each rule, the program
  1. enumerates the spatially periodic configurations B of spatial period p <= PMAX that satisfy the
     necessary condition for a stuck right edge of reach1.py (same routine);
  2. for 1 <= b <= BMAX and 0 < |a| <= b builds the window graph of Fix(sigma^a F^b) (vertices: the
     words of length 2b; edges: the words w of length 2b+1 with F^b(w) = w_{b+a}), pruned to the
     vertices on bi-infinite paths;
  3. for every kept B that lies in Fix(sigma^a F^b), that is, all of whose windows of length 2b+1
     are edges, checks that every window of B of length 2b has exactly one outgoing edge in the
     pruned graph. Then the windows of B form a cycle that no path leaves, and B is the only point
     of Fix(sigma^a F^b) that agrees with B on a left half-line.
It prints, per rule, the number of kept backgrounds, the number of triples (B, a, b) with B in
Fix(sigma^a F^b), and the number of those in which some window of B has a second outgoing edge.
Needs NumPy only.

Usage: python reach1_unique.py RULE [RULE ...]      (defaults PMAX = 6, BMAX = 4)
       python reach1_unique.py classes              the 44 rules of the 11 classes of the proposition
       python reach1_unique.py controls             rules 110, 14, 93, 25, 87, 159
"""
import sys
import numpy as np

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
    """list of (p, tuple of cells) for the kept periodic configurations (as in reach1.py)"""
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
        for _ in range(2 * p + 2):
            mask &= mask[Q]; Q = Q[Q]
        for s in np.nonzero(mask)[0]: out.append((p, tuple(int(v) for v in C[s])))
    return out

def fix_graph(rule, a, b):
    """edge set (as a boolean vector over the words of length 2b+1) and out-degrees of the pruned graph"""
    L = 2 * b + 1; idx = np.arange(4 ** L, dtype=np.int64); W = cells_of(idx, L); V = W
    for _ in range(b): V = step(rule, V)
    e = idx[V[:, 0] == W[:, b + a]]
    n = 4 ** (2 * b)
    while True:
        src = e >> 2; dst = e & (n - 1)
        hs = np.zeros(n, dtype=bool); hp = np.zeros(n, dtype=bool); hs[src] = True; hp[dst] = True
        alive = hs & hp
        k = alive[src] & alive[dst]
        if k.all(): break
        e = e[k]
    is_edge = np.zeros(4 ** L, dtype=bool); is_edge[e] = True
    return is_edge, np.bincount(e >> 2, minlength=n)

def run(rule):
    B = backgrounds(rule); tested = 0; bad = []
    for b in range(1, BMAX + 1):
        for a in list(range(-b, 0)) + list(range(1, b + 1)):
            is_edge, outdeg = fix_graph(rule, a, b)
            for (p, cfg) in B:
                edges = {sum(cfg[(j + i) % p] << (2 * (2 * b - i)) for i in range(2 * b + 1)) for j in range(p)}
                if not all(is_edge[w] for w in edges): continue
                tested += 1
                if any(outdeg[w >> 2] != 1 for w in edges): bad.append((a, b, p, ''.join(map(str, cfg))))
    print(f"rule {rule}: {len(B)} backgrounds kept, {tested} triples (B,a,b) with B in Fix(sigma^a F^b), "
          f"{len(bad)} with a second continuation" + (f"; first: {bad[:5]}" if bad else ""), flush=True)
    return len(bad)

if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["classes"]:
        rules = sum((cls(r) for r in (28, 29, 41, 43, 47, 57, 61, 111, 143, 155, 157)), [])
    elif args == ["controls"]:
        rules = [110, 14, 93, 25, 87, 159]
    else:
        rules = [int(x) for x in args]
    total = sum(run(r) for r in rules)
    print("total number of triples with a second continuation:", total)
