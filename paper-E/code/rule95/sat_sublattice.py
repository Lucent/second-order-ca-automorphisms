"""Block maps h of radius R commuting with F_95 whose output bits read one sublattice each ("sublattice type"):
     type P (par=0): y'_i = TY[b], b_m = y_{i+m} (m even), x_{i+m} (m odd);  x'_i = TX[b], b_m = x_{i+m} (m even), y_{i+m} (m odd)
     type S (par=1): the same with even and odd exchanged;        index = sum_m b_m << (m+R), |m| <= R.
   Unknowns: the 2*2^(2R+1) table entries.  The condition hF = Fh at cell 0 is written out for all values of the
   2R+3 bits (cells -R-1..R+1) on which each output bit of both sides depends.
   mode 'all':  enumerate all solutions, classify.
   mode 'aut':  enumerate the restrictions to Phi^+ (and to Phi^-) that occur among the solutions, keep the injective
                ones, and enumerate all solutions with an injective restriction on both; an automorphism of sublattice
                type and radius <= R is among them.
   usage: sat_sublattice.py R all|aut [SOLVER [RULE]]     (needs the package python-sat; default solver cadical153)
   The result is written to sat_R<R>_<mode>.pkl in the current directory.
   RULE (default 95) may be any elementary rule f that does not read its middle argument; the condition is then
   TY[r0] = TX[r1] xor f(TY[r2], ., TY[r3]).  For RULE other than 95 only the mode 'all' is available: all solutions
   are enumerated and classified as sigma^a F^b, constant or other, each solution is rechecked against hF = Fh on
   random configurations, and the result is written to sat_rule<RULE>_R<R>_all.pkl."""
import sys, time, itertools, pickle, numpy as np
from enum_golden import gwords, injective

RULE = 95
FT = None; G2 = None
def set_rule(rule):
    """f = elementary rule `rule`, which must not read its middle argument; G2[a][c] = f(a, ., c)"""
    global RULE, FT, G2
    f = lambda a, b, c: (rule >> (4 * a + 2 * b + c)) & 1
    assert all(f(a, 0, c) == f(a, 1, c) for a in (0, 1) for c in (0, 1)), "the rule reads its middle argument"
    RULE = rule
    FT = np.array([[[f(a, b, c) for c in (0, 1)] for b in (0, 1)] for a in (0, 1)], dtype=np.int64)
    G2 = [[f(a, 0, c) for c in (0, 1)] for a in (0, 1)]
set_rule(95)

class Enc:
    def __init__(self, R, par):
        self.R, self.par, self.N = R, par, 1 << (2 * R + 1)
    def vx(self, w): return 1 + w
    def vy(self, w): return 1 + self.N + w
    def idx(self, X, Y, i, out):
        """index of the window of output bit `out` ('x' or 'y') at cell i; X, Y: dict cell -> array/int"""
        R = self.R; w = 0
        for m in range(-R, R + 1):
            same = (m - self.par) % 2 == 0
            lay = (Y if same else X) if out == 'y' else (X if same else Y)
            w = w | (lay[i + m] << (m + R))
        return w
    def class_cells(self, cl, lo, hi, other):
        """arrays over all assignments of the bits of class cl (class 0: y_k k even, x_k k odd) in cells lo..hi"""
        n = hi - lo + 1; a = np.arange(1 << n, dtype=np.int64); X = {}; Y = {}
        for j, k in enumerate(range(lo, hi + 1)):
            b = (a >> j) & 1; o = np.full(1 << n, other, dtype=np.int64)
            if (k % 2 == 0) == (cl == 0): Y[k], X[k] = b, o
            else: X[k], Y[k] = b, o
        return X, Y
    def constraints(self):
        R = self.R; eqs = set(); ands = set(); found = {'x': 0, 'y': 0}
        for cl in (0, 1):
            res = []
            for other in (0, 1):
                X, Y = self.class_cells(cl, -R - 1, R + 1, other)
                FX = {k: Y[k] for k in range(-R, R + 1)}
                FY = {k: X[k] ^ FT[Y[k - 1], Y[k], Y[k + 1]] for k in range(-R, R + 1)}
                res.append((self.idx(FX, FY, 0, 'y'), self.idx(X, Y, 0, 'x'), self.idx(X, Y, -1, 'y'), self.idx(X, Y, 1, 'y'),
                            self.idx(FX, FY, 0, 'x'), self.idx(X, Y, 0, 'y')))
            ydep = all((res[0][j] == res[1][j]).all() for j in range(4)); xdep = all((res[0][j] == res[1][j]).all() for j in (4, 5))
            assert ydep != xdep
            r = res[0]
            if ydep:    # TY[r0] = TX[r1] xor f(TY[r2], ., TY[r3]);  for f = 95: TX[r1] xor 1 xor TY[r2] TY[r3]
                found['y'] += 1
                ands |= set(zip((self.vy(r[0])).tolist(), (self.vx(r[1])).tolist(), (self.vy(r[2])).tolist(), (self.vy(r[3])).tolist()))
            else:       # TX[r4] = TY[r5]
                found['x'] += 1
                eqs |= set(zip((self.vx(r[4])).tolist(), (self.vy(r[5])).tolist()))
        assert found == {'x': 1, 'y': 1}
        cnf = []
        for a, b in eqs:
            if a != b: cnf += [[a, -b], [-a, b]]
        for t in ands:
            vs = sorted(set(t))
            for vals in itertools.product((0, 1), repeat=len(vs)):
                d = dict(zip(vs, vals))
                if d[t[0]] != (d[t[1]] ^ G2[d[t[2]]][d[t[3]]]):
                    cnf.append([-v if d[v] else v for v in vs])
        return cnf, len(eqs), len(ands)
    # ----- restrictions to Phi^+ / Phi^- as pairs of block maps of G (window R+1): lists of (G-word, variable)
    def restriction(self, sign):
        R = self.R; W = gwords(R + 1); out = {}
        for ob in ('x', 'y'):
            hit = None
            for cls in (0, 1):
                qs = [q for q in range(-R - 1, R + 1) if q % 2 == cls]      # R+1 values, ascending
                rows = []
                for other in (1,):
                    lst = []
                    for a in W:
                        d = {q: 1 for q in range(-R - 2, R + 2)}
                        for j, q in enumerate(qs): d[q] = 1 - ((a >> (R - j)) & 1)
                        if sign > 0: X = {i: d[-1 - i] for i in range(-R, R + 1)}; Y = {i: d[-i] for i in range(-R, R + 1)}
                        else: X = {i: d[i - 1] for i in range(-R, R + 1)}; Y = {i: d[i] for i in range(-R, R + 1)}
                        w = self.idx(X, Y, 0, ob); lst.append((a, self.vx(w) if ob == 'x' else self.vy(w)))
                    rows.append(lst)
                # dependent class: the one whose variation changes the index
                if len(set(v for _, v in rows[0])) > 1: assert hit is None; hit = rows[0]
            assert hit is not None; out[ob] = hit
        return out

def named(R, enc):
    """tables (TX,TY) of sigma^a F^b, |a|+|b| <= R, of the type of enc; dict key -> name"""
    names = {}
    n = 2 * R + 1
    for a in range(-R, R + 1):
        for b in range(-R + abs(a), R - abs(a) + 1):
            if (a + b) % 2 != enc.par: continue
            tabs = []
            for ob in ('x', 'y'):
                # the class read by this output bit: find it via the index function
                T = np.zeros(enc.N, dtype=np.uint8); done = False
                for cl in (0, 1):
                    X, Y = enc.class_cells(cl, -R, R, 0); w = enc.idx(X, Y, 0, ob)
                    if len(np.unique(w)) != enc.N: continue
                    X2, Y2 = enc.class_cells(cl, -R, R, 0)
                    lo, hi = -R, R
                    for _ in range(abs(b)):
                        if b > 0: NX = {k: Y2[k] for k in range(lo + 1, hi)}; NY = {k: X2[k] ^ FT[Y2[k - 1], Y2[k], Y2[k + 1]] for k in range(lo + 1, hi)}
                        else: NY = {k: X2[k] for k in range(lo + 1, hi)}; NX = {k: Y2[k] ^ FT[X2[k - 1], X2[k], X2[k + 1]] for k in range(lo + 1, hi)}
                        X2, Y2 = NX, NY; lo += 1; hi -= 1
                    val = (X2 if ob == 'x' else Y2)[a]          # (sigma^a c)_0 = c_a
                    T[w] = val; done = True
                assert done; tabs.append(T)
            names[tabs[0].tobytes() + tabs[1].tobytes()] = (a, b)
    return names

def model_key(model, enc):
    m = np.array(model[:2 * enc.N]) > 0
    return m.astype(np.uint8).tobytes()

def apply_table(R, par, TX, TY, x, y):
    """the block map with tables TX, TY of type par on a configuration of period n (arrays x, y)"""
    n = len(x); wx = np.zeros(n, dtype=np.int64); wy = np.zeros(n, dtype=np.int64)
    for m in range(-R, R + 1):
        same = (m - par) % 2 == 0
        xs, ys = np.roll(x, -m), np.roll(y, -m)
        wy |= (ys if same else xs) << (m + R); wx |= (xs if same else ys) << (m + R)
    return TX[wx].astype(np.int64), TY[wy].astype(np.int64)

def step(x, y):
    return y, x ^ FT[np.roll(y, 1), y, np.roll(y, -1)]

def run_rule(R, par, solver_name="cadical153", log=sys.stdout):
    """mode 'all' for a rule other than 95: enumerate, classify, recheck"""
    from pysat.solvers import Solver
    t0 = time.time(); enc = Enc(R, par); cnf, neq, nand = enc.constraints(); nv = 2 * enc.N
    print("rule %d R=%d type %s: %d variables, %d equalities, %d constraints of four entries, %d clauses [%.1fs]" % (
        RULE, R, "PS"[par], nv, neq, nand, len(cnf), time.time() - t0), file=log, flush=True)
    names = named(R, enc)
    s = Solver(name=solver_name, bootstrap_with=cnf); sols = []
    while s.solve():
        m = s.get_model(); m = m + [-(i) for i in range(len(m) + 1, nv + 1)]
        sols.append(m[:nv]); s.add_clause([-l for l in m[:nv]])
    s.delete()
    rng = np.random.default_rng(1); n = 4 * R + 9; bad = 0
    cnt = {"named": 0, "const": 0, "other": 0}; rows = []; seen = set()
    for m in sols:
        key = model_key(m, enc)
        if key in names: cnt["named"] += 1; lab = "s^%d F^%d" % names[key]; seen.add(names[key])
        elif len(set(key)) == 1: cnt["const"] += 1; lab = "const"
        else: cnt["other"] += 1; lab = "other"
        t = np.frombuffer(key, dtype=np.uint8); TX, TY = t[:enc.N], t[enc.N:]
        for _ in range(20):
            x = rng.integers(0, 2, n, dtype=np.int64); y = rng.integers(0, 2, n, dtype=np.int64)
            a = apply_table(R, par, TX, TY, *step(x, y)); b = step(*apply_table(R, par, TX, TY, x, y))
            bad += not ((a[0] == b[0]).all() and (a[1] == b[1]).all())
        rows.append((lab, key))
    print("  all solutions: %d; %s; maps sigma^a F^b of this type: %d, found: %d; failures of the recheck: %d [%.1fs]" % (
        len(sols), cnt, len(names), len(seen), bad, time.time() - t0), file=log, flush=True)
    return {"R": R, "par": par, "rule": RULE, "rows": rows, "counts": cnt}

def run(R, par, mode, solver_name="cadical153", log=sys.stdout):
    from pysat.solvers import Solver
    t0 = time.time(); enc = Enc(R, par); cnf, neq, nand = enc.constraints()
    nv = 2 * enc.N
    print("R=%d type %s: %d variables, %d equalities, %d and-constraints, %d clauses [%.1fs]" % (R, "PS"[par], nv, neq, nand, len(cnf), time.time() - t0), file=log, flush=True)
    names = named(R, enc)
    allv = list(range(1, nv + 1))
    res = {"R": R, "par": par}
    def restr_info(model, rs):
        val = lambda v: 1 if model[v - 1] > 0 else 0
        info = []
        for ob in ('x', 'y'):
            f = {a: 1 - val(v) for a, v in rs[ob]}
            info.append(tuple(f[a] for a in gwords(R + 1)))
        return tuple(info)
    def is_inj(info):
        W = gwords(R + 1)
        return all(injective(R + 1, dict(zip(W, t))) for t in info)
    rsP, rsM = enc.restriction(+1), enc.restriction(-1)
    if mode == "all":
        s = Solver(name=solver_name, bootstrap_with=cnf); sols = []
        while s.solve():
            m = s.get_model(); m = m + [-(i) for i in range(len(m) + 1, nv + 1)]
            sols.append(m[:nv]); s.add_clause([-l for l in m[:nv]])
        s.delete()
        cnt = {"named": 0, "const": 0, "other": 0, "other_injP": 0, "other_injM": 0, "other_inj_both": 0}
        rows = []
        for m in sols:
            key = model_key(m, enc); iP, iM = restr_info(m, rsP), restr_info(m, rsM)
            if key in names: cnt["named"] += 1; lab = "s^%d F^%d" % names[key]
            elif len(set(key)) == 1: cnt["const"] += 1; lab = "const"
            else:
                a, b = is_inj(iP), is_inj(iM)
                cnt["other"] += 1; cnt["other_injP"] += a; cnt["other_injM"] += b; cnt["other_inj_both"] += a and b; lab = "other"
            rows.append((lab, key, iP, iM))
        print("  all solutions: %d; %s; named expected %d; distinct (Phi+,Phi-) restrictions %d, distinct Phi+ %d, distinct Phi- %d [%.1fs]" % (
            len(sols), cnt, len(names), len(set((r[2], r[3]) for r in rows)), len(set(r[2] for r in rows)), len(set(r[3] for r in rows)), time.time() - t0), file=log, flush=True)
        res["rows"] = rows; res["counts"] = cnt
    else:
        inj = {}
        for sign, rs in ((+1, rsP), (-1, rsM)):
            pv = sorted(set(v for ob in rs for _, v in rs[ob]))
            s = Solver(name=solver_name, bootstrap_with=cnf); seen = []
            while s.solve():
                m = s.get_model(); m = m + [-(i) for i in range(len(m) + 1, nv + 1)]
                seen.append((restr_info(m, rs), [m[v - 1] for v in pv])); s.add_clause([-m[v - 1] for v in pv])
                if len(seen) % 2000 == 0: print("    ... %d restrictions [%.0fs]" % (len(seen), time.time() - t0), file=log, flush=True)
            s.delete()
            good = [x for x in seen if is_inj(x[0])]
            inj[sign] = good
            print("  restrictions to Phi^%s among the solutions: %d (projection on %d variables); injective: %d [%.1fs]" % ("+" if sign > 0 else "-", len(seen), len(pv), len(good), time.time() - t0), file=log, flush=True)
            res["nrestr%+d" % sign] = len(seen); res["restr%+d" % sign] = [x[0] for x in seen]
        s = Solver(name=solver_name, bootstrap_with=cnf); found = []
        for a in inj[+1]:
            for b in inj[-1]:
                ass = a[1] + b[1]; k = 0
                while s.solve(assumptions=ass):
                    m = s.get_model(); m = m + [-(i) for i in range(len(m) + 1, nv + 1)]
                    key = model_key(m, enc); found.append((names.get(key), key, a[0], b[0])); k += 1
                    s.add_clause([-l for l in m[:nv]])
        s.delete()
        nn = sum(1 for x in found if x[0] is not None)
        print("  solutions with injective restrictions to Phi^+ and Phi^-: %d, of which sigma^a F^b: %d (expected named of this type: %d), others: %d [%.1fs]" % (
            len(found), nn, len(names), len(found) - nn, time.time() - t0), file=log, flush=True)
        res["found"] = found
    return res

if __name__ == "__main__":
    R = int(sys.argv[1]); mode = sys.argv[2]; solver = sys.argv[3] if len(sys.argv) > 3 else "cadical153"
    rule = int(sys.argv[4]) if len(sys.argv) > 4 else 95
    out = {}
    if rule == 95:
        for par in (0, 1): out[par] = run(R, par, mode, solver)
        pickle.dump(out, open("sat_R%d_%s.pkl" % (R, mode), "wb"))
    else:
        assert mode == "all", "for a rule other than 95 only the mode 'all' is available"
        set_rule(rule)
        for par in (0, 1): out[par] = run_rule(R, par, solver)
        tot = {k: out[0]["counts"][k] + out[1]["counts"][k] for k in out[0]["counts"]}
        print("rule %d R=%d: both types together %s; 2R^2+2R+1 = %d" % (rule, R, tot, 2 * R * R + 2 * R + 1), flush=True)
        pickle.dump(out, open("sat_rule%d_R%d_all.pkl" % (rule, R), "wb"))
