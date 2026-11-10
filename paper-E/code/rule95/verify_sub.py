"""For the solutions of sublattice type found by sat_sublattice.py (mode 'all', radius R): each non-constant solution is
   hatH_{phi,psi} sigma^n with particle-conserving block maps phi, psi: (1) restrictions to Phi^+, Phi^- of the form
   of hatH_{phi,psi} sigma^n (exact), phi, psi conserving (exact); (2) the table generated from the definition of
   hatH_{phi,psi} sigma^n by simulation equals the solution in all entries.   usage: verify_sub.py R sat_R<R>_<mode>.pkl"""
import sys, pickle, numpy as np
from core import *
from enum_golden import gwords, injective
from gen7 import build
from census_cmp import conserving, same_map, shifted, minimal

def class_maps(TX, TY, par, R):
    M = {}
    for sign in (+1, -1):
        for q0, ob, t in ((0, 'y', TY), (-1, 'x', TX)):
            cls = (q0 % 2) ^ par
            qs = [q for q in range(-R - 1, R + 1) if q % 2 == cls]; js = [(q - cls) // 2 for q in qs]
            f = {}
            for a in gwords(len(qs)):
                d = {q: 1 for q in range(-R - 2, R + 2)}
                for j, q in enumerate(qs): d[q] = 1 - ((a >> (len(qs) - 1 - j)) & 1)
                if sign > 0: cx = {i: d[-1 - i] for i in range(-R, R + 1)}; cy = {i: d[-i] for i in range(-R, R + 1)}
                else: cx = {i: d[i - 1] for i in range(-R, R + 1)}; cy = {i: d[i] for i in range(-R, R + 1)}
                w = 0
                for m in range(-R, R + 1):
                    same = (m - par) % 2 == 0
                    lay = (cy if same else cx) if ob == 'y' else (cx if same else cy)
                    w |= lay[m] << (m + R)
                f[a] = 1 - int(t[w])
            jout = (q0 - (q0 % 2)) // 2
            M[(sign, q0)] = GMap(len(qs), js[0] - jout, f)
    if par == 0:
        phi, psi = M[(1, 0)], M[(-1, 0)]; form = same_map(M[(1, -1)], phi) and same_map(M[(-1, -1)], psi)
    else:
        phi, psi = M[(1, -1)], M[(-1, 0)]; form = same_map(M[(1, 0)], shifted(phi, -1)) and same_map(M[(-1, -1)], shifted(psi, +1))
    return phi, psi, form

if __name__ == "__main__":
    R = int(sys.argv[1]); src = sys.argv[2]
    d = pickle.load(open(src, "rb")); N = 1 << (2 * R + 1)
    for par in (0, 1):
        rows = d[par]["rows"] if "rows" in d[par] else [("other", x[1], None, None) for x in d[par]["found"] if x[0] is None]
        st = {"n": 0, "const": 0, "form": 0, "conserving": 0, "inj": 0, "table_equal": 0}
        for k, row in enumerate(rows):
            key = row[1]; a = np.frombuffer(key, dtype=np.uint8); TX, TY = a[:N], a[N:]
            if len(set(key)) == 1: st["const"] += 1; continue
            st["n"] += 1
            phi, psi, form = class_maps(TX, TY, par, R)
            cons = conserving(phi.w, phi.tab) and conserving(psi.w, psi.tab)
            st["form"] += form; st["conserving"] += cons; st["inj"] += injective(phi.w, phi.tab) and injective(psi.w, psi.tab)
            if form and cons:
                BX, BY, inc = build(phi, psi, nsig=par, R=R, par_y=par, log=False, seed=k)
                st["table_equal"] += bool(inc == 0 and (BX == TX).all() and (BY == TY).all())
            if "found" in d[par]: print("   ", "PS"[par], "phi:", minimal(phi)[:14], "inj", injective(phi.w, phi.tab), "| psi:", minimal(psi)[:14], "inj", injective(psi.w, psi.tab), flush=True)
        print("R=%d type %s:" % (R, "PS"[par]), st, flush=True)
