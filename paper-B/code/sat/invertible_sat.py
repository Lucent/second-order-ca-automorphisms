"""Invertible part of the centralizer by SAT: same commutation encoding as centralizer.py plus Hedlund's
balance condition (an invertible, hence surjective, block map takes every output symbol on exactly
4^{2R} of the 4^{2R+1} table entries), encoded as cardinality constraints. Every balanced solution is then
tested for injectivity with the pair-graph test. Output: the exact set of invertible block maps of radius R
commuting with F (the balance condition is necessary, the pair-graph test is exact)."""
import sys, time, json, argparse
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType
from centralizer import *

def run(rule, R, maxsol, native=False):
    f, F_loc, Finv_loc = make_rule(rule)
    n = 2 * R + 1; N = 4 ** n; nvar = 2 * N
    names = named_maps(R, F_loc, Finv_loc)
    with Solver(name="minicard" if native else "cadical153") as s:
        t0 = time.time()
        nv, nc = encode(R, f, F_loc, s)
        top = nvar
        # indicator variables I_{e,sym} <-> (X_e, Y_e) == sym, then sum_e I_{e,sym} == 4^{2R} for each sym
        for sym in range(4):
            xb, yb = sym >> 1, sym & 1
            inds = []
            for e in range(N):
                top += 1; I = top; X, Y = 2 * e + 1, 2 * e + 2
                lx = X if xb else -X; ly = Y if yb else -Y
                s.add_clause([-I, lx]); s.add_clause([-I, ly]); s.add_clause([I, -lx, -ly])
                inds.append(I)
            if native:   # minicard: sum(I) <= k and sum(not I) <= n-k  <=>  sum(I) == k
                s.add_atmost(inds, 4 ** (2 * R)); s.add_atmost([-i for i in inds], N - 4 ** (2 * R))
            else:
                cnf = CardEnc.equals(lits=inds, bound=4 ** (2 * R), top_id=top, encoding=EncType.seqcounter)
                top = max(top, cnf.nv)
                for cl in cnf.clauses: s.add_clause(cl)
        t_enc = time.time() - t0; t0 = time.time()
        sols = []; capped = False
        while s.solve():
            model = s.get_model()
            val = [1 if lit > 0 else 0 for lit in model[:nvar]]
            sols.append([(val[2 * e] << 1) | val[2 * e + 1] for e in range(N)])
            s.add_clause([-lit for lit in model[:nvar]])
            if len(sols) % 100 == 0: print(f"    ... {len(sols)} balanced solutions [{time.time()-t0:.0f}s]", file=sys.stderr, flush=True)
            if len(sols) >= maxsol: capped = True; break
        t_enum = time.time() - t0
    print(f"\n=== rule {rule}R radius {R}: balanced commuting maps: {len(sols)}{' (CAPPED)' if capped else ''}; "
          f"encode {t_enc:.1f}s, enumeration {t_enum:.1f}s", flush=True)
    rows = []; ninv = 0; nother_inv = 0
    for tab in sols:
        g = BlockMap(R, tab)
        assert check_commutes(g, F_loc, R)
        name = names.get(g.key())
        inj = True if name is not None else injective_pairgraph_np(g)
        label = f"sigma^{name[0]} F^{name[1]}" if name else "OTHER"
        rows.append({"label": label, "invertible": inj, "layer": layer_structure(g)[0], "affine": is_affine(g), "table": g.table_str()})
        ninv += inj; nother_inv += (inj and not name)
    print(f"  invertible: {ninv}; of which sigma^a F^b: {ninv - nother_inv} (expected {len(names)}); OTHER invertible: {nother_inv}")
    for r in sorted(rows, key=lambda r: (r["label"] == "OTHER", r["label"])):
        if r["invertible"] or len(rows) <= 40:
            print(f"    {r['label']:14s} inv={r['invertible']!s:5s} {r['layer']:15s} affine={r['affine']}  table={r['table'][:80]}{'...' if len(r['table'])>80 else ''}")
    return {"rule": rule, "R": R, "n_balanced": len(sols), "capped": capped, "n_invertible": ninv, "n_other_invertible": nother_inv,
            "t_encode": t_enc, "t_enum": t_enum, "solutions": rows}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rules", default="54"); ap.add_argument("--R", default="2")
    ap.add_argument("--max", type=int, default=100000); ap.add_argument("--out", default=None); ap.add_argument("--native", action="store_true")
    a = ap.parse_args(); res = []
    for R in map(int, a.R.split(",")):
        for rule in map(int, a.rules.split(",")):
            res.append(run(rule, R, a.max, a.native))
            if a.out: json.dump(res, open(a.out, "w"), indent=1)
