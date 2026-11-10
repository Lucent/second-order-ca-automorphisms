"""Exact comparison of the census tables of 95R (all block maps of radius R commuting with F_95; radius <= 3:
   census_R3.json; radius 4: a text file, one table per line) with the maps hatH_{phi,psi} sigma^n, phi, psi particle-conserving block maps of G.
   For every census solution: (1) recheck hF = Fh on all 4^(2R+3) windows; (2) each output bit reads one sublattice
   (exact, all 4^(2R+1) entries); (3) the restrictions to Phi^+ and Phi^- are phi, psi on both classes, in the form of
   hatH_{phi,psi} sigma^n with n in {0,1} (exact, all words); phi, psi conserve particles (coboundary criterion) and
   map G into G; injectivity of phi, psi (exact); (4) the table generated from the definition of hatH_{phi,psi} sigma^n
   by simulation equals the census table in all entries.
   usage: census_cmp.py R [FILE]
   R <= 3: FILE (default census_R3.json) is a JSON list with one entry per radius, each with the fields "R" and
   "solutions", a list of objects with the fields "label" and "table".  R >= 4: FILE (default census_R<R>.txt) has one
   table per line.  A table is a string of 4^(2R+1) digits 0..3: the digit at position e is the image cell 2x'+y' of
   the window whose cells c_{-R}, ..., c_R, each written 2x+y, are the base-4 digits of e, c_{-R} most significant."""
import sys, os, json, numpy as np
from core import *
from enum_golden import gwords, injective
from gen7 import build

def conserving(w, f):
    """f: dict G-word(len w) -> 0/1; exact coboundary test, and image in G"""
    if w == 1: return f[0] == 0 and f.get(1, 1) == 1
    mask = (1 << (w - 1)) - 1; lam = {0: 0}; stack = [0]
    E = gwords(w)
    adj = {}
    for e in E:
        s, t = e >> 1, e & mask; wt = f[e] - ((e >> (w - 1)) & 1)        # f(e) - e_0 = lam(s) - lam(t)
        adj.setdefault(s, []).append((t, -wt)); adj.setdefault(t, []).append((s, wt))
    while stack:
        v = stack.pop()
        for u, d in adj[v]:
            if u not in lam: lam[u] = lam[v] + d; stack.append(u)
            elif lam[u] != lam[v] + d: return False
    for e in E:
        for b in (0, 1):
            e2 = ((e & mask) << 1) | b
            if e2 in f and f[e] and f[e2]: return False
    return True

def extend(g, lo, hi):
    """table of the GMap g on the window [lo,hi] (relative offsets), as dict over G-words"""
    assert lo <= g.a and g.a + g.w - 1 <= hi
    W = hi - lo + 1; out = {}
    for m in gwords(W):
        sub = (m >> (hi - (g.a + g.w - 1))) & ((1 << g.w) - 1); out[m] = g.tab[sub]
    return out
def same_map(g1, g2):
    lo = min(g1.a, g2.a); hi = max(g1.a + g1.w, g2.a + g2.w) - 1
    return extend(g1, lo, hi) == extend(g2, lo, hi)
def shifted(g, d): return GMap(g.w, g.a + d, g.tab)

def analyse(tab, R, verbose=False):
    n = 2 * R + 1; N = 4 ** n
    T = np.frombuffer(tab.encode(), dtype=np.uint8) - 48
    assert len(T) == N
    e = np.arange(N, dtype=np.int64)
    sym = {m: (e >> (2 * (R - m))) & 3 for m in range(-R, R + 1)}       # cell m; first cell most significant
    X = {m: sym[m] >> 1 for m in sym}; Y = {m: sym[m] & 1 for m in sym}
    res = {}
    for par in (0, 1):
        ok = True; tabs = []
        for ob in ('x', 'y'):
            w = np.zeros(N, dtype=np.int64)
            for m in range(-R, R + 1):
                same = (m - par) % 2 == 0
                lay = (Y if same else X) if ob == 'y' else (X if same else Y)
                w |= lay[m] << (m + R)
            val = (T >> 1) if ob == 'x' else (T & 1)
            t = np.full(1 << n, 255, dtype=np.uint8); t[w] = val
            if not (t[w] == val).all(): ok = False; break
            tabs.append(t)
        if ok: res[par] = tabs
    if not res: return {"type": None}
    const = len(set(T.tolist())) == 1
    par = 0 if 0 in res else 1
    if const: return {"type": "const", "value": int(T[0])}
    TX, TY = res[par]
    # restrictions: four class maps
    M = {}
    for sign in (+1, -1):
        for q0, ob, t in ((0, 'y', TY), (-1, 'x', TX)):
            # the cell 0 of the image is (pi'(-1), pi'(0)) resp. (xi'(-1), xi'(0))
            for cls in (0, 1):
                qs = [q for q in range(-R - 1, R + 1) if q % 2 == cls]; js = [(q - cls) // 2 for q in qs]
                f = {}
                for a in gwords(len(qs)):
                    vals = set()
                    for other in (0, 1):               # the other class carries no particle
                        if other: continue
                        d = {q: 1 for q in range(-R - 2, R + 2)}
                        for j, q in enumerate(qs): d[q] = 1 - ((a >> (len(qs) - 1 - j)) & 1)
                        if sign > 0: cx = {i: d[-1 - i] for i in range(-R, R + 1)}; cy = {i: d[-i] for i in range(-R, R + 1)}
                        else: cx = {i: d[i - 1] for i in range(-R, R + 1)}; cy = {i: d[i] for i in range(-R, R + 1)}
                        w = 0
                        for m in range(-R, R + 1):
                            same = (m - par) % 2 == 0
                            lay = (cy if same else cx) if ob == 'y' else (cx if same else cy)
                            w |= lay[m] << (m + R)
                        vals.add(1 - int(t[w]))
                    f[a] = vals.pop()
                if len(set(f.values())) > 1 or cls == ((q0 % 2) ^ par):
                    pass
                # the class that is read: output class is q0 mod 2; type P reads the same class, type S the other
                if cls == ((q0 % 2) ^ par):
                    jout = (q0 - (q0 % 2)) // 2
                    M[(sign, q0)] = GMap(len(qs), js[0] - jout, f)
    out = {"type": "PS"[par], "TX": TX, "TY": TY}
    if par == 0:
        phi, psi = M[(1, 0)], M[(-1, 0)]
        form = same_map(M[(1, -1)], phi) and same_map(M[(-1, -1)], psi)
    else:
        phi, psi = M[(1, -1)], M[(-1, 0)]
        form = same_map(M[(1, 0)], shifted(phi, -1)) and same_map(M[(-1, -1)], shifted(psi, +1))
    out.update(phi=phi, psi=psi, form=form, cons=(conserving(phi.w, phi.tab), conserving(psi.w, psi.tab)),
               inj=(injective(phi.w, phi.tab), injective(psi.w, psi.tab)))
    return out

def commutes(tab, R):
    T = np.frombuffer(tab.encode(), dtype=np.uint8).astype(np.int64) - 48
    n = 2 * R + 1; m = n + 2; idx = np.arange(4 ** m, dtype=np.int64)
    c = [(idx >> (2 * (m - 1 - i))) & 3 for i in range(m)]
    def Floc(s0, s1, s2):
        x, y = s1 >> 1, s1 & 1; a, cc = s0 & 1, s2 & 1
        return (y << 1) | (x ^ ((95 >> (4 * a + 2 * y + cc)) & 1))
    Fw = [Floc(c[i], c[i + 1], c[i + 2]) for i in range(n)]
    lhs = T[sum(Fw[i] << (2 * (n - 1 - i)) for i in range(n))]
    g = [T[sum(c[j + i] << (2 * (n - 1 - i)) for i in range(n))] for j in range(3)]
    return bool((lhs == Floc(g[0], g[1], g[2])).all())

def minimal(g):
    """minimal window description of a GMap (string) for the log"""
    return "w=%d a=%d tab=%s" % (g.w, g.a, "".join(str(g.tab[e]) for e in gwords(g.w)))

if __name__ == "__main__":
    R = int(sys.argv[1])
    HERE = os.path.dirname(os.path.abspath(__file__))
    if R <= 3:
        d = json.load(open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "census_R3.json")))
        sols = [(s["label"], s["table"]) for x in d if x["R"] == R for s in x["solutions"]]
    else:
        sols = [("?", l.strip()) for l in open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "census_R%d.txt" % R)) if l.strip()]
    stats = {"total": len(sols), "commute": 0, "const": 0, "sublattice": 0, "form": 0, "conserving": 0, "both_injective": 0, "noninjective": 0, "table_equal": 0, "P": 0, "S": 0}
    pairs = set()
    for k, (label, tab) in enumerate(sols):
        stats["commute"] += commutes(tab, R)
        a = analyse(tab, R)
        if a["type"] == "const": stats["const"] += 1; continue
        if a["type"] is None: print("solution", k, label, "NOT of sublattice type"); continue
        stats["sublattice"] += 1; stats[a["type"]] += 1
        stats["form"] += a["form"]; stats["conserving"] += all(a["cons"])
        if all(a["inj"]): stats["both_injective"] += 1
        else: stats["noninjective"] += 1
        par = 0 if a["type"] == "P" else 1
        TX, TY, inc = build(a["phi"], a["psi"], nsig=par, R=R, par_y=par, log=False, seed=k)
        eq = inc == 0 and (TX == a["TX"]).all() and (TY == a["TY"]).all()
        stats["table_equal"] += bool(eq)
        pairs.add((minimal(a["phi"]), minimal(a["psi"]), par))
        if not (a["form"] and all(a["cons"]) and eq): print("solution", k, label, "FAILS:", a["form"], a["cons"], eq, inc)
        if R == 3 and not all(a["inj"]): print("  ", label, a["type"], "phi:", minimal(a["phi"]), "inj", a["inj"][0], "| psi:", minimal(a["psi"]), "inj", a["inj"][1])
    print("R=%d:" % R, stats, "distinct (phi,psi,n):", len(pairs))
