"""Exact small-radius centralizer of a second-order reversible CA among ALL block maps on the full
4-symbol alphabet, by SAT.

Automaton. Alphabet A = {0,1}^2, cell s = (x, y) encoded as the integer s = 2x + y (y = current row,
x = previous row). For an elementary rule f (radius 1),
    F(c)_i = ( y_i ,  x_i XOR f(y_{i-1}, y_i, y_{i+1}) ).
F is reversible with inverse F^{-1}(c)_i = ( y_i XOR f(x_{i-1}, x_i, x_{i+1}) , x_i ).

Unknown. A block map g : A^{2R+1} -> A given by its lookup table: entry e in [0, 4^(2R+1)) is the
(2R+1)-word read as a base-4 integer (leftmost symbol most significant); g[e] = (X_e, Y_e) are two
boolean unknowns, SAT variables 2e+1 (X) and 2e+2 (Y).

Commutation. g o F and F o g are both block maps of radius R+1; they are equal iff their local rules
agree on every window w in A^{2R+3}:
    (g o F)(w)_center = g[ F(w) ]            -- F(w) is a concrete (2R+1)-word, so this is ONE table entry eL
    (F o g)(w)_center = F( g[e1], g[e2], g[e3] )  -- e1,e2,e3 = the three (2R+1)-subwindows of w
                     = ( x = Y_{e2} ,  y = X_{e2} XOR f(Y_{e1}, Y_{e2}, Y_{e3}) ).
So each window gives the two constraints
    X_{eL} = Y_{e2}                              (2 binary clauses)
    Y_{eL} XOR X_{e2} = f(Y_{e1}, Y_{e2}, Y_{e3})   (<= 16 clauses, by truth table over the <= 5 variables)
Entries may coincide (e.g. eL = e2); the truth-table encoding handles that by working over the set of
distinct variables. Total: 4^(2R+3) windows, 2*4^(2R+1) variables, <= 18 clauses per window.
All solutions are enumerated with blocking clauses (cap --max).
"""
import sys, time, json, argparse, itertools
from pysat.solvers import Solver

# ------------------------------------------------------------------ the automaton
def make_rule(rule):
    f = lambda a, b, c: (rule >> ((a << 2) | (b << 1) | c)) & 1
    def F_loc(s0, s1, s2):
        x, y = s1 >> 1, s1 & 1
        return (y << 1) | (x ^ f(s0 & 1, y, s2 & 1))
    def Finv_loc(s0, s1, s2):
        x, y = s1 >> 1, s1 & 1
        return ((y ^ f(s0 >> 1, x, s2 >> 1)) << 1) | x
    return f, F_loc, Finv_loc

def apply_local(loc, word, r):
    """apply a radius-r local rule (function of 2r+1 symbols) to a word; output has length len-2r."""
    k = 2 * r + 1
    return tuple(loc(*word[i:i + k]) for i in range(len(word) - k + 1))

def word_index(word):
    e = 0
    for s in word: e = (e << 2) | s
    return e

def index_word(e, n):
    return tuple((e >> (2 * (n - 1 - i))) & 3 for i in range(n))

# ------------------------------------------------------------------ block maps as tables (radius R, 4^(2R+1) entries)
class BlockMap:
    """table[e] = output symbol for the (2R+1)-word with index e (centred neighbourhood [-R, R])."""
    def __init__(self, R, table):
        self.R = R; self.table = tuple(table); self.n = 2 * R + 1
    def __call__(self, word):            # apply to a word, output shorter by 2R
        k = self.n
        return tuple(self.table[word_index(word[i:i + k])] for i in range(len(word) - k + 1))
    def key(self): return self.table
    @staticmethod
    def from_local(R, loc, r):
        """local rule of radius r <= R embedded as a radius-R table."""
        n = 2 * R + 1; tab = []
        for e in range(4 ** n):
            w = index_word(e, n)
            tab.append(loc(*w[R - r:R + r + 1]))
        return BlockMap(R, tab)
    @staticmethod
    def shift(R, k):
        """sigma^k: g(w) = w_{centre + k}, i.e. the output at i is c_{i+k}."""
        n = 2 * R + 1
        return BlockMap(R, [index_word(e, n)[R + k] for e in range(4 ** n)])
    @staticmethod
    def const(R, s):
        return BlockMap(R, [s] * 4 ** (2 * R + 1))
    def compose(self, other, R):
        """self o other, as a radius-R table (requires self.R + other.R <= R)."""
        n = 2 * R + 1; r1, r2 = self.R, other.R; tab = []
        for e in range(4 ** n):
            w = index_word(e, n)
            inner = w[R - r1 - r2:R + r1 + r2 + 1]
            mid = other(inner)                        # length 2 r1 + 1
            tab.append(self.table[word_index(mid)] if mid else None)
        return BlockMap(R, tab)
    def table_str(self):
        return "".join(str(s) for s in self.table)

def named_maps(R, F_loc, Finv_loc):
    """all sigma^a F^b with |a| + |b| <= R, as radius-R tables, keyed by table."""
    out = {}
    for b in range(-R, R + 1):
        rb = abs(b)
        if b == 0: Fb = BlockMap.from_local(rb, lambda *w: w[0], 0)
        else:
            loc = F_loc if b > 0 else Finv_loc
            Fb = BlockMap.from_local(1, loc, 1)
            for _ in range(rb - 1): Fb = BlockMap.from_local(1, loc, 1).compose(Fb, Fb.R + 1)
        for a in range(-(R - rb), R - rb + 1):
            g = BlockMap.shift(R - rb, a).compose(Fb, R) if R - rb >= 0 else None
            out[g.key()] = (a, b)
    return out

# ------------------------------------------------------------------ SAT encoding
def encode(R, f, F_loc, solver):
    n = 2 * R + 1; m = 2 * R + 3
    X = lambda e: 2 * e + 1
    Y = lambda e: 2 * e + 2
    nclauses = 0
    for ew in range(4 ** m):
        w = index_word(ew, m)
        eL = word_index(apply_local(F_loc, w, 1))
        e1, e2, e3 = word_index(w[0:n]), word_index(w[1:n + 1]), word_index(w[2:n + 2])
        # new x-bit: X_eL = Y_e2
        solver.add_clause([-X(eL), Y(e2)]); solver.add_clause([X(eL), -Y(e2)]); nclauses += 2
        # new y-bit: Y_eL = X_e2 xor f(Y_e1, Y_e2, Y_e3)
        roles = [Y(eL), X(e2), Y(e1), Y(e2), Y(e3)]
        vs = sorted(set(roles))
        for bits in itertools.product((0, 1), repeat=len(vs)):
            val = dict(zip(vs, bits))
            a = [val[v] for v in roles]
            if (a[0] ^ a[1]) != f(a[2], a[3], a[4]):
                solver.add_clause([-v if val[v] else v for v in vs]); nclauses += 1
    return 2 * 4 ** n, nclauses

def enumerate_solutions(R, f, F_loc, maxsol, solver_name="cadical153"):
    n = 2 * R + 1; nvar = 2 * 4 ** n
    with Solver(name=solver_name) as s:
        t0 = time.time()
        nv, nc = encode(R, f, F_loc, s)
        t_enc = time.time() - t0
        sols = []; t0 = time.time(); capped = False
        while s.solve():
            model = s.get_model()
            val = [1 if lit > 0 else 0 for lit in model[:nvar]]
            table = [(val[2 * e] << 1) | val[2 * e + 1] for e in range(4 ** n)]
            sols.append(table)
            s.add_clause([-lit for lit in model[:nvar]])
            if len(sols) % 1000 == 0:
                print(f"    ... {len(sols)} solutions after {time.time()-t0:.0f}s", file=sys.stderr, flush=True)
            if len(sols) >= maxsol: capped = True; break
        t_enum = time.time() - t0
    return sols, capped, nv, nc, t_enc, t_enum

# ------------------------------------------------------------------ verification and classification
def check_commutes(g, F_loc, R):
    """direct re-check of g o F = F o g on all windows (independent of the SAT encoding)."""
    m = 2 * R + 3
    for ew in range(4 ** m):
        w = index_word(ew, m)
        if g(apply_local(F_loc, w, 1)) != apply_local(F_loc, g(w), 1): return False
    return True

def injective_pairgraph(g):
    """Injectivity of g on A^Z by the pair (de Bruijn product) graph: nodes (u,v) in A^{2R} x A^{2R},
    edge (u,v)->(u',v') on symbols (a,b) iff g(u a) = g(v b). g is non-injective iff a bi-infinite path
    passes through an off-diagonal node; bi-infinite paths live in the subgraph left after iteratively
    deleting nodes with no in- or no out-edge."""
    R = g.R; L = 2 * R; N = 4 ** L; mask = N - 1; tab = g.table
    succ = [None] * (N * N)
    for u in range(N):
        for v in range(N):
            lst = []
            for a in range(4):
                gu = tab[(u << 2) | a]
                for b in range(4):
                    if tab[(v << 2) | b] == gu:
                        lst.append((((u << 2) | a) & mask) * N + (((v << 2) | b) & mask))
            succ[u * N + v] = lst
    alive = [True] * (N * N)
    indeg = [0] * (N * N)
    for i in range(N * N):
        for j in succ[i]: indeg[j] += 1
    outdeg = [len(l) for l in succ]
    stack = [i for i in range(N * N) if indeg[i] == 0 or outdeg[i] == 0]
    pred = [[] for _ in range(N * N)]
    for i in range(N * N):
        for j in succ[i]: pred[j].append(i)
    while stack:
        i = stack.pop()
        if not alive[i]: continue
        alive[i] = False
        for j in succ[i]:
            indeg[j] -= 1
            if alive[j] and indeg[j] == 0: stack.append(j)
        for j in pred[i]:
            outdeg[j] -= 1
            if alive[j] and outdeg[j] == 0: stack.append(j)
    for u in range(N):
        for v in range(N):
            if u != v and alive[u * N + v]: return False
    return True

def injective_pairgraph_np(g):
    """same pair-graph test, vectorised with numpy (nodes = pairs of 2R-words; 16 candidate edges each)."""
    import numpy as np
    R = g.R; L = 2 * R; N = 4 ** L; mask = N - 1
    tab = np.array(g.table, dtype=np.int8)
    u = (np.arange(N * N, dtype=np.int32) // N).astype(np.int32); v = (np.arange(N * N, dtype=np.int32) % N).astype(np.int32)
    a = np.arange(4, dtype=np.int32); b = np.arange(4, dtype=np.int32)
    gu = tab[(u[:, None] << 2) | a[None, :]]                 # (N², 4)
    gv = tab[(v[:, None] << 2) | b[None, :]]
    ok_out = (gu[:, :, None] == gv[:, None, :]).reshape(N * N, 16)
    su = ((u[:, None] << 2) | a[None, :]) & mask; sv = ((v[:, None] << 2) | b[None, :]) & mask
    succ = (su[:, :, None] * N + sv[:, None, :]).reshape(N * N, 16)
    # predecessors of (u,v): (c u>>2, d v>>2), edge symbols (u&3, v&3)
    c = np.arange(4, dtype=np.int32); d = np.arange(4, dtype=np.int32)
    pu = (c[None, :] << (2 * (L - 1))) | (u[:, None] >> 2); pv = (d[None, :] << (2 * (L - 1))) | (v[:, None] >> 2)
    gpu = tab[(pu[:, :, None] << 2) | (u & 3)[:, None, None]]; gpv = tab[(pv[:, None, :] << 2) | (v & 3)[:, None, None]]
    ok_in = (gpu == gpv).reshape(N * N, 16)
    pred = (pu[:, :, None] * N + pv[:, None, :]).reshape(N * N, 16)
    alive = np.ones(N * N, dtype=bool)
    while True:
        has_out = (ok_out & alive[succ]).any(axis=1)
        has_in = (ok_in & alive[pred]).any(axis=1)
        new = alive & has_out & has_in
        if (new == alive).all(): break
        alive = new
    offdiag = u != v
    return not (alive & offdiag).any()

def injective_periodic(g, maxperiod):
    """cross-check: injective on spatially periodic configurations of period <= maxperiod."""
    R = g.R
    for p in range(1, maxperiod + 1):
        seen = {}
        for e in range(4 ** p):
            c = index_word(e, p)
            ext = (c * (2 * R // p + 3))[:p + 2 * R]
            # rotate so that the window is centred: output at position i uses c[i-R..i+R]
            ext = tuple(c[(i - R) % p] for i in range(p + 2 * R))
            img = g(ext)
            if img in seen and seen[img] != c: return False
            seen[img] = c
    return True

def layer_structure(g):
    """Returns ('layerwise-same', G) if g(x,y) = (G x, G y) for one binary block map G;
    ('layer-separable', (G1,G2)) if g(x,y) = (G1 x, G2 y); else ('joint', None).
    Also reports 'cross' if g(x,y) = (G1 y, G2 x)."""
    R = g.R; n = 2 * R + 1
    def project(bit_out, bit_in):
        """does output bit depend only on input layer bit_in? if so return its truth table."""
        tabl = {}
        for e in range(4 ** n):
            w = index_word(e, n)
            key = tuple((s >> bit_in) & 1 for s in w)
            val = (g.table[e] >> bit_out) & 1
            if tabl.setdefault(key, val) != val: return None
        return tuple(tabl[k] for k in itertools.product((0, 1), repeat=n))
    Gx = project(1, 1); Gy = project(0, 0)          # X_out from x-layer, Y_out from y-layer
    if Gx is not None and Gy is not None:
        return ("layerwise-same" if Gx == Gy else "layer-separable", (Gx, Gy))
    Gxy = project(1, 0); Gyx = project(0, 1)
    if Gxy is not None and Gyx is not None: return ("cross", (Gxy, Gyx))
    return ("joint", None)

def is_affine(g):
    """is g affine over F_2 (as a map (F_2^2)^{2R+1} -> F_2^2)? check g(w)+g(0) linear."""
    n = 2 * R_of(g); tab = g.table; z = tab[0]
    N = 4 ** (2 * g.R + 1)
    # linearity: g(u xor v) = g(u) xor g(v) xor g(0) for basis vectors suffices if we check all pairs with basis
    basis = [1 << i for i in range(2 * (2 * g.R + 1))]
    for e in range(N):
        for bvec in basis:
            if tab[e ^ bvec] ^ z != (tab[e] ^ z) ^ (tab[bvec] ^ z): return False
    return True
def R_of(g): return g.R

# ------------------------------------------------------------------ main
def run(rule, R, maxsol, do_inj=True):
    f, F_loc, Finv_loc = make_rule(rule)
    # sanity: F o F^{-1} = id, F^{-1} o F = id (as radius-2 tables)
    Fm = BlockMap.from_local(1, F_loc, 1); Fi = BlockMap.from_local(1, Finv_loc, 1)
    idm = BlockMap.shift(2, 0)
    assert Fm.compose(Fi, 2).key() == idm.key() and Fi.compose(Fm, 2).key() == idm.key(), "F not reversible?!"
    names = named_maps(R, F_loc, Finv_loc)
    sols, capped, nv, nc, t_enc, t_enum = enumerate_solutions(R, f, F_loc, maxsol)
    print(f"\n=== rule {rule}R, radius R={R}: {nv} variables, {nc} clauses, encode {t_enc:.1f}s, "
          f"enumeration {t_enum:.1f}s, {len(sols)} solutions{' (CAPPED)' if capped else ''}", flush=True)
    rows = []
    t0 = time.time()
    counts = {"named": 0, "const": 0, "other": 0, "invertible": 0, "layerwise-same": 0,
              "layer-separable": 0, "cross": 0, "joint": 0, "affine": 0}
    for table in sols:
        g = BlockMap(R, table)
        assert check_commutes(g, F_loc, R), "SAT solution does not commute: encoding bug"
        name = names.get(g.key())
        if name is not None:
            a, b = name; label = f"sigma^{a} F^{b}"; counts["named"] += 1
        elif len(set(table)) == 1:
            label = f"const ({table[0] >> 1},{table[0] & 1})"; counts["const"] += 1
        else:
            label = "OTHER"; counts["other"] += 1
        inv = None
        if do_inj:
            if name is not None: inv = True
            elif len(set(table)) == 1: inv = False
            else:
                inv = injective_pairgraph_np(g)
                if R <= 1: assert inv == injective_pairgraph(g), "numpy/python pair-graph mismatch"
                if R <= 2 and inv: assert injective_periodic(g, 2 * R + 4), "pair-graph says injective but periodic test fails"
            if inv: counts["invertible"] += 1
        kind, _ = layer_structure(g); counts[kind] += 1
        aff = is_affine(g); counts["affine"] += int(aff)
        rows.append({"label": label, "invertible": inv, "layer": kind, "affine": aff, "table": g.table_str()})
    t_cls = time.time() - t0
    # all named maps must have been found
    missing = [v for k, v in names.items() if k not in {tuple(int(c) for c in r["table"]) for r in rows}]
    assert not missing or capped, f"named maps missing from solutions: {missing}"
    print(f"  classification {t_cls:.1f}s: {counts}; expected named sigma^a F^b with |a|+|b|<=R: {len(names)}; "
          f"all named found: {not missing}")
    rows.sort(key=lambda r: (r["label"] == "OTHER", r["label"]))
    if len(rows) <= 64:
        for r in rows:
            print(f"    {r['label']:18s} inv={str(r['invertible']):5s} {r['layer']:16s} affine={r['affine']}  table={r['table']}")
    else:
        for r in rows:
            if r["label"] != "OTHER" or counts["other"] <= 64:
                print(f"    {r['label']:18s} inv={str(r['invertible']):5s} {r['layer']:16s} affine={r['affine']}" +
                      (f"  table={r['table']}" if R <= 1 else ""))
    return {"rule": rule, "R": R, "nvars": nv, "nclauses": nc, "t_encode": t_enc, "t_enum": t_enum,
            "t_classify": t_cls, "nsol": len(sols), "capped": capped, "counts": counts, "n_named_expected": len(names),
            "solutions": rows}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", default="110,30,90,54")
    ap.add_argument("--R", default="1,2")
    ap.add_argument("--max", type=int, default=100000)
    ap.add_argument("--out", default=None)
    ap.add_argument("--noinj", action="store_true")
    args = ap.parse_args()
    results = []
    for R in [int(x) for x in args.R.split(",")]:
        for rule in [int(x) for x in args.rules.split(",")]:
            results.append(run(rule, R, args.max, not args.noinj))
            if args.out:
                with open(args.out, "w") as fh: json.dump(results, fh, indent=1)
