"""The 36 maps hatH_{phi,psi} sigma of gen36.py: exact verification of the generated radius-7 tables (tables36.npz).
   (A) hF = Fh on all 2*2^17 class windows; (B) restrictions to Phi^+ and Phi^- (all words) are those of
   hatH_{phi,psi} sigma.  A block map commuting with F is determined by its restrictions to Phi^+ and Phi^-, so the
   table map is then hatH_{phi,psi} sigma, which so has radius <= 7.
   (C) exact neighbourhood of each.  (D) if sat_others_R7_S.npz (SAT search, written by sat_post.py) is present: the
   two sets of tables coincide.   usage: verify36.py [tables36.npz [sat_others_R7_S.npz]]"""
import numpy as np, sys, os
from core import *
from sat_sublattice import Enc, FT
from verify_sub import class_maps
from census_cmp import same_map
from gen36 import PHI, PSI
R = 7; enc = Enc(R, 1); N = enc.N
HERE = os.path.dirname(os.path.abspath(__file__))
f36 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "tables36.npz")
fsat = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "sat_others_R7_S.npz")
gen = np.load(f36)
def commutes(TX, TY):
    bad = 0
    for cl in (0, 1):
        res = []
        for other in (0, 1):
            X, Y = enc.class_cells(cl, -R - 1, R + 1, other)
            FX = {k: Y[k] for k in range(-R, R + 1)}; FY = {k: X[k] ^ FT[Y[k - 1], Y[k], Y[k + 1]] for k in range(-R, R + 1)}
            gFy = TY[enc.idx(FX, FY, 0, 'y')]; gFx = TX[enc.idx(FX, FY, 0, 'x')]
            gy = {k: TY[enc.idx(X, Y, k, 'y')] for k in (-1, 0, 1)}; gx0 = TX[enc.idx(X, Y, 0, 'x')]
            Fgx = gy[0]; Fgy = gx0 ^ FT[gy[-1], gy[0], gy[1]]
            res.append(((gFy != Fgy), (gFx != Fgx)))
        # the bit that depends on this class only: mismatch pattern must be identical for both fillings; count all mismatches
        bad += int(res[0][0].sum() + res[0][1].sum() + res[1][0].sum() + res[1][1].sum())
    return bad
tabs = {}; allok = True
for na, a in PHI:
    for nb, b in PSI:
        key = "%s | %s" % (na, nb)
        if key not in gen.files: continue
        TX = gen[key][:N].astype(np.int64); TY = gen[key][N:].astype(np.int64)
        bad = commutes(TX, TY)
        phi, psi, form = class_maps(TX, TY, 1, R)
        ok = bad == 0 and form and same_map(phi, a) and same_map(psi, b)
        idx = np.arange(N)
        depy = [m for m in range(-R, R + 1) if (TY != TY[idx ^ (1 << (m + R))]).any()]; depx = [m for m in range(-R, R + 1) if (TX != TX[idx ^ (1 << (m + R))]).any()]
        nb_lo, nb_hi = min(depy + depx), max(depy + depx)
        tabs[(na, nb)] = TX.astype(np.uint8).tobytes() + TY.astype(np.uint8).tobytes()
        allok &= ok
        print("%-8s | %-12s | commutation mismatches %d, restrictions as predicted: %s, neighbourhood [%d,%d]" % (na, nb, bad, form and same_map(phi, a) and same_map(psi, b), nb_lo, nb_hi), flush=True)
print("tables verified:", len(tabs), "all ok:", allok, "; pairwise distinct:", len(set(tabs.values())) == len(tabs))
if os.path.exists(fsat):
    S = set(row.tobytes() for row in np.load(fsat)["tables"])
    print("SAT search: %d automorphisms of sublattice type S and radius <= 7 outside <sigma,F>; equal to the set of the %d tables above: %s" % (len(S), len(tabs), S == set(tabs.values())))
