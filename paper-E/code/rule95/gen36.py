"""Which hatH_{phi,psi} sigma, with phi in {sigma^-3..sigma^4, phi_* sigma, phibar_*} and psi in {sigma^-4..sigma^3, phi_*,
   phibar_* sigma^-1} (the injective maps allowed by the window lemma at radius 7), not both shifts, have radius <= 7?
   For each pair the table of radius 7 is filled from the definition by simulation, and every sample must be consistent
   with the entries already filled.  The consistent tables are written to tables36.npz (key "phi | psi", value the 2^15
   entries of TX followed by those of TY); verify36.py checks them exactly.   usage: gen36.py [outfile]"""
import numpy as np, sys, os
from core import *
from gen7 import build
phi1 = gm_swap('0001001', '0000101'); phi1m = gm_swap('1001000', '1010000')
sh = lambda g, d: GMap(g.w, g.a + d, g.tab)        # g o sigma^d
PHI = [("s^%d" % k, gm_shift(k)) for k in range(-3, 5)] + [("phi* s", sh(phi1, 1)), ("phibar*", phi1m)]
PSI = [("s^%d" % k, gm_shift(k)) for k in range(-4, 4)] + [("phi*", phi1), ("phibar* s^-1", sh(phi1m, -1))]
if __name__ == '__main__':
    res = {}; tabs = {}
    for na, a in PHI:
        for nb, b in PSI:
            if na.startswith("s^") and nb.startswith("s^"): continue
            try:
                TX, TY, inc = build(a, b, nsig=1, R=7, par_y=1, log=False)
            except AssertionError as e:
                inc = -1
            res[(na, nb)] = inc
            print(na, "|", nb, "| inconsistent samples:", inc, flush=True)
            if inc == 0: tabs["%s | %s" % (na, nb)] = np.concatenate((TX, TY)).astype(np.uint8)
    print("consistent with radius 7:", sorted(k for k, v in res.items() if v == 0))
    np.savez_compressed(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables36.npz"), **tabs)
