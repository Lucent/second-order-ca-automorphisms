"""Balanced commuting maps (invertible-part candidates) by SAT, memory-lean version of invertible_sat.py:
commutation clauses (centralizer.encode) + Hedlund balance as native minicard at-most constraints; blocking
clauses use one representative SAT variable per union-find class (X_eL = Y_e2 equalities), and every solution's
table is appended to raw_balanced_<rule>_R<R>.txt as soon as it is found."""
import sys, time, json, argparse
import numpy as np
from pysat.solvers import Solver
from centralizer import *
from propsearch import build

def run(rule, R, maxsol):
    f, F_loc, Finv_loc = make_rule(rule)
    n = 2 * R + 1; N = 4 ** n; nvar = 2 * N
    S = build(rule, R); cls = S["cls"]
    reps = {}
    for bit in range(nvar): reps.setdefault(int(cls[bit]), bit + 1)       # SAT var = bit index + 1
    rep_vars = sorted(reps.values())
    print(f"rule {rule}R R={R}: {nvar} table bits, {len(rep_vars)} classes", flush=True)
    out = open(f"raw_balanced_{rule}_R{R}.txt", "w")
    with Solver(name="minicard") as s:
        t0 = time.time()
        nv, nc = encode(R, f, F_loc, s)
        top = nvar
        for sym in range(4):
            xb, yb = sym >> 1, sym & 1; inds = []
            for e in range(N):
                top += 1; I = top; X, Y = 2 * e + 1, 2 * e + 2
                lx = X if xb else -X; ly = Y if yb else -Y
                s.add_clause([-I, lx]); s.add_clause([-I, ly]); s.add_clause([I, -lx, -ly]); inds.append(I)
            s.add_atmost(inds, 4 ** (2 * R)); s.add_atmost([-i for i in inds], N - 4 ** (2 * R))
        t_enc = time.time() - t0; t0 = time.time(); nsol = 0
        while s.solve():
            model = s.get_model()
            val = [1 if lit > 0 else 0 for lit in model[:nvar]]
            out.write("".join(str((val[2 * e] << 1) | val[2 * e + 1]) for e in range(N)) + "\n"); out.flush()
            nsol += 1
            s.add_clause([-model[v - 1] for v in rep_vars])
            if nsol % 500 == 0: print(f"    ... {nsol} balanced solutions [{time.time()-t0:.0f}s]", flush=True)
            if nsol >= maxsol: print("CAPPED"); break
        t_enum = time.time() - t0
    out.close()
    print(f"=== rule {rule}R radius {R}: balanced commuting maps: {nsol}; encode {t_enc:.1f}s, enumeration {t_enum:.1f}s", flush=True)
    json.dump({"rule": rule, "R": R, "n_balanced": nsol, "capped": nsol >= maxsol, "t_encode": t_enc, "t_enum": t_enum,
               "tables_file": f"raw_balanced_{rule}_R{R}.txt"}, open(f"raw_balanced_{rule}_R{R}.json", "w"))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--rule", type=int, default=54); ap.add_argument("--R", type=int, default=3)
    ap.add_argument("--max", type=int, default=200000)
    a = ap.parse_args(); run(a.rule, a.R, a.max)
