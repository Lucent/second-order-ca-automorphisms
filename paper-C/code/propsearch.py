"""Memory-light exact enumeration of the radius-R centralizer by constraint propagation + backtracking (numpy).
Same constraints as centralizer.py:  for each window w in A^{2R+3}:
    X_eL = Y_e2                                   -> merged by union-find into one variable class
    Y_eL xor X_e2 xor f(Y_e1, Y_e2, Y_e3) = 0     -> one 5-slot constraint over classes
Search: DFS over unknown classes; after each decision, propagate: every constraint with exactly one unknown class
is evaluated for both values of that class (forced value / conflict). A complete assignment is verified against all
constraints and decoded into a lookup table (checked again with check_commutes_np-style evaluation).
Memory: O(#windows) int32, i.e. ~100 MB at R=4 (versus ~12 GB for the 75M-clause CNF)."""
import sys, time, json, argparse
import numpy as np
from centralizer import make_rule, named_maps, BlockMap, layer_structure, is_affine, injective_pairgraph_np

def build(rule, R):
    f, F_loc, Finv_loc = make_rule(rule)
    n = 2 * R + 1; m = 2 * R + 3; N = 4 ** n; W = 4 ** m
    idx = np.arange(W, dtype=np.int64)
    cells = [((idx >> (2 * (m - 1 - i))) & 3).astype(np.int8) for i in range(m)]
    Ft = np.array([[[F_loc(a, b, c) for c in range(4)] for b in range(4)] for a in range(4)], dtype=np.int8)
    Fw = [Ft[cells[i], cells[i + 1], cells[i + 2]] for i in range(m - 2)]
    eL = sum(Fw[i].astype(np.int64) << (2 * (n - 1 - i)) for i in range(n))
    e1 = (idx >> 4) & (N - 1); e2 = (idx >> 2) & (N - 1); e3 = idx & (N - 1)
    del cells, Fw
    # union-find over variables: X_e -> 2e, Y_e -> 2e+1 ; merge X_eL ~ Y_e2
    nv = 2 * N
    parent = np.arange(nv, dtype=np.int64)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    pairs = np.unique(np.stack([2 * eL, 2 * e2 + 1], 1), axis=0)
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb: parent[max(ra, rb)] = min(ra, rb)
    # full path compression vectorised
    for _ in range(64):
        p2 = parent[parent]
        if (p2 == parent).all(): break
        parent = p2
    roots, cls = np.unique(parent, return_inverse=True)
    ncls = len(roots)
    cons = np.stack([cls[2 * eL + 1], cls[2 * e2], cls[2 * e1 + 1], cls[2 * e2 + 1], cls[2 * e3 + 1]], 1).astype(np.int32)
    cons = np.unique(cons, axis=0)
    ftab = np.array([f(a, b, c) for a in range(2) for b in range(2) for c in range(2)], dtype=np.int8)
    # watch lists (CSR): class -> constraint ids
    flat_cls = cons.reshape(-1); flat_con = np.repeat(np.arange(len(cons), dtype=np.int32), 5)
    order = np.argsort(flat_cls, kind="stable")
    wl_con = flat_con[order]; wl_ptr = np.searchsorted(flat_cls[order], np.arange(ncls + 1))
    return dict(f=f, F_loc=F_loc, Finv_loc=Finv_loc, n=n, N=N, ncls=ncls, cls=cls, cons=cons, ftab=ftab,
                wl_con=wl_con, wl_ptr=wl_ptr)

def evaluate(v, ftab):
    return (v[:, 0] ^ v[:, 1] ^ ftab[(v[:, 2].astype(np.int32) << 2) | (v[:, 3].astype(np.int32) << 1) | v[:, 4].astype(np.int32)]) == 0

def propagate(S, val, active):
    cons, ftab = S["cons"], S["ftab"]
    while len(active):
        c = cons[active]; v = val[c]
        unk = v < 0; nunk = unk.sum(1)
        full = nunk == 0
        if full.any() and not evaluate(v[full], ftab).all(): return False
        cand = nunk > 0
        if not cand.any(): return True
        cc = c[cand]; vv = v[cand]; uu = unk[cand]
        big = np.where(uu, cc, np.iinfo(np.int32).max).min(1); small = np.where(uu, cc, -1).max(1)
        single = big == small
        if not single.any(): return True
        cc, vv, uu, u = cc[single], vv[single], uu[single], big[single]
        v0 = np.where(uu, 0, vv).astype(np.int8); v1 = np.where(uu, 1, vv).astype(np.int8)
        s0 = evaluate(v0, ftab); s1 = evaluate(v1, ftab)
        if (~s0 & ~s1).any(): return False
        force = s0 ^ s1
        if not force.any(): return True
        fu = u[force]; fv = s1[force].astype(np.int8)
        # contradictions among forced values
        order = np.argsort(fu, kind="stable"); fu, fv = fu[order], fv[order]
        first = np.r_[True, fu[1:] != fu[:-1]]
        if len(fu) > 1:
            # same class forced to different values?
            grp_min = np.minimum.reduceat(fv, np.flatnonzero(first)); grp_max = np.maximum.reduceat(fv, np.flatnonzero(first))
            if (grp_min != grp_max).any(): return False
        fu, fv = fu[first], fv[first]
        val[fu] = fv
        # new active: constraints watching the changed classes
        starts, ends = S["wl_ptr"][fu], S["wl_ptr"][fu + 1]
        active = np.unique(np.concatenate([S["wl_con"][s:e] for s, e in zip(starts, ends)]))
    return True

def search(S, maxsol, verbose=True):
    ncls = S["ncls"]; sols = []; t0 = time.time(); nodes = 0
    val = -np.ones(ncls, dtype=np.int8)
    if not propagate(S, val, np.arange(len(S["cons"]), dtype=np.int32)): return sols, nodes
    stack = [val]
    while stack:
        val = stack.pop(); nodes += 1
        unknown = np.flatnonzero(val < 0)
        if len(unknown) == 0:
            if evaluate(val[S["cons"]], S["ftab"]).all(): sols.append(val.copy())
            if verbose and len(sols) % 10 == 0: print(f"    ... {len(sols)} solutions, {nodes} nodes [{time.time()-t0:.0f}s]", flush=True)
            if len(sols) >= maxsol: break
            continue
        # branch on the unknown class with the most constraints
        deg = S["wl_ptr"][unknown + 1] - S["wl_ptr"][unknown]
        u = unknown[np.argmax(deg)]
        for b in (1, 0):
            v2 = val.copy(); v2[u] = b
            act = S["wl_con"][S["wl_ptr"][u]:S["wl_ptr"][u + 1]]
            if propagate(S, v2, act): stack.append(v2)
    return sols, nodes

def decode(S, val, R):
    cls = S["cls"]; N = S["N"]
    X = val[cls[2 * np.arange(N)]]; Y = val[cls[2 * np.arange(N) + 1]]
    return [int((X[e] << 1) | Y[e]) for e in range(N)]

def run(rule, R, maxsol, do_inj):
    t0 = time.time(); S = build(rule, R); t_build = time.time() - t0
    print(f"\n=== rule {rule}R radius {R} (propagation search): {len(S['cons'])} constraints over {S['ncls']} variable classes "
          f"(from {2*S['N']} table bits); build {t_build:.1f}s", flush=True)
    t0 = time.time(); sols, nodes = search(S, maxsol); t_search = time.time() - t0
    print(f"  search: {len(sols)} solutions, {nodes} search nodes, {t_search:.1f}s", flush=True)
    names = named_maps(R, S["F_loc"], S["Finv_loc"])
    from classify_r4 import check_commutes_np_R
    counts = {"named": 0, "const": 0, "other": 0, "invertible": 0, "layerwise-same": 0, "joint": 0, "layer-separable": 0, "cross": 0}
    rows = []
    for v in sols:
        tab = decode(S, v, R)
        assert check_commutes_np_R(tab, S["F_loc"], R), "decoded solution fails independent commutation check"
        g = BlockMap(R, tab); name = names.get(g.key())
        if name: label = f"sigma^{name[0]} F^{name[1]}"; counts["named"] += 1; inv = True
        elif len(set(tab)) == 1: label = f"const ({tab[0]>>1},{tab[0]&1})"; counts["const"] += 1; inv = False
        else:
            label = "OTHER"; counts["other"] += 1
            inv = injective_pairgraph_np(g) if (do_inj and R <= 3) else None
        counts["invertible"] += bool(inv)
        kind = layer_structure(g)[0]; counts[kind] += 1
        rows.append({"label": label, "invertible": inv, "layer": kind})
    missing = [v for k, v in names.items() if k not in {tuple(decode(S, v, R)) for v in sols}]
    print(f"  classification: {counts}; expected named {len(names)}; all named found: {not missing}")
    for r in sorted(rows, key=lambda r: (r["label"] == "OTHER", r["label"])): print(f"    {r['label']:16s} inv={r['invertible']} {r['layer']}")
    return {"rule": rule, "R": R, "ncons": int(len(S["cons"])), "ncls": int(S["ncls"]), "t_build": t_build, "t_search": t_search,
            "nodes": nodes, "nsol": len(sols), "counts": counts, "solutions": rows}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rules", default="110"); ap.add_argument("--R", default="1,2,3")
    ap.add_argument("--max", type=int, default=100000); ap.add_argument("--out", default=None); ap.add_argument("--noinj", action="store_true")
    a = ap.parse_args(); res = []
    for R in map(int, a.R.split(",")):
        for rule in map(int, a.rules.split(",")):
            res.append(run(rule, R, a.max, not a.noinj))
            if a.out: json.dump(res, open(a.out, "w"), indent=1)
