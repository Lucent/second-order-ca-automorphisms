"""Count of the particle-conserving block maps of G with window w <= 7 by linear algebra:
   conservation is imposed as the linear equations  sum_e n_e(z) f(e) = |z|  over all cyclic G-words z of length
   <= #vertices (every simple cycle of the de Bruijn graph), the affine solution space is parametrised by free
   coordinates, and all 0/1 values of the free coordinates are tried; the result is compared with the depth-first
   search of enum_golden.py.  Also identifies the two non-shift injective maps of window 8 with phi_* (exchange of
   0001001 and 0000101) and phibar_* (exchange of 1001000 and 1010000).   usage: enum_golden_linear.py WMAX"""
import sys, numpy as np
from fractions import Fraction
from enum_golden import gwords, enumerate_maps, injective, table_str
from core import gm_swap

def cyc_words(n):   # cyclic G-words of length n (as tuples), all rotations included (harmless)
    return [m for m in range(1 << n) if not (m & (m >> 1)) and not ((m & 1) and (m >> (n - 1)) & 1)] if n > 1 else [0]

def brute(w):
    E = gwords(w); ei = {e: i for i, e in enumerate(E)}; nV = len(gwords(w - 1))
    rows = set()
    for n in range(1, nV + 1):
        for m in cyc_words(n):
            bits = [(m >> (n - 1 - i)) & 1 for i in range(n)]
            cnt = [0] * len(E)
            for i in range(n):
                e = 0
                for k in range(w): e = (e << 1) | bits[(i + k) % n]
                cnt[ei[e]] += 1
            rows.add(tuple(cnt) + (sum(bits),))
    # Gaussian elimination over Q
    basis = {}                                         # pivot col -> row (list of Fractions)
    for r in sorted(rows):
        r = [Fraction(x) for x in r]
        for pc, br in basis.items():
            if r[pc]: 
                c = r[pc]; r = [a - c * b for a, b in zip(r, br)]
        piv = next((i for i in range(len(E)) if r[i]), None)
        if piv is None:
            assert r[-1] == 0; continue
        c = r[piv]; r = [a / c for a in r]
        for pc in list(basis):
            br = basis[pc]
            if br[piv]:
                c = br[piv]; basis[pc] = [a - c * b for a, b in zip(br, r)]
        basis[piv] = r
    free = [i for i in range(len(E)) if i not in basis]
    nf = len(free)
    den = 1
    from math import lcm
    for r in basis.values():
        for a in r: den = lcm(den, a.denominator)
    idx = np.arange(1 << nf, dtype=np.int64)
    fv = np.zeros((1 << nf, len(E)), dtype=np.int64)
    for k, c in enumerate(free): fv[:, c] = (idx >> k) & 1
    ok = np.ones(1 << nf, dtype=bool)
    for pc, r in basis.items():
        val = np.full(1 << nf, int(r[-1] * den), dtype=np.int64)
        for c in free:
            if r[c]: val -= int(r[c] * den) * fv[:, c]
        ok &= (val == 0) | (val == den)
        fv[:, pc] = val // den
    fv = fv[ok]
    # image in G: f(e) f(e') = 0 for consecutive edges
    mask = (1 << (w - 1)) - 1
    for e in E:
        for b in (0, 1):
            e2 = ((e & mask) << 1) | b
            if e2 in ei: fv = fv[(fv[:, ei[e]] & fv[:, ei[e2]]) == 0]
    return len(rows), len(basis), nf, set("".join(map(str, row)) for row in fv)

for w in range(2, int(sys.argv[1]) + 1):
    nrows, rank, nf, S = brute(w)
    sols, _ = enumerate_maps(w)
    S2 = set(table_str(w, f) for f in sols)
    print("w=%d: %d cycle equations, rank %d, free %d; brute count %d; DFS count %d; equal sets: %s" % (w, nrows, rank, nf, len(S), len(S2), S == S2), flush=True)

# window 8: the non-shift injective maps
sols, _ = enumerate_maps(8); E = gwords(8)
inj = [f for f in sols if injective(8, f)]
phi1 = gm_swap('0001001', '0000101'); phi1m = gm_swap('1001000', '1010000')
print("phi_* window", phi1.w, "anchor", phi1.a, "; phibar_* window", phi1m.w, "anchor", phi1m.a)
t1 = "".join(str(phi1.tab[e]) for e in E); t2 = "".join(str(phi1m.tab[e]) for e in E)
names = []
for f in inj:
    s = table_str(8, f)
    names.append("sigma^4 phi_*" if s == t1 else "sigma^3 phibar_*" if s == t2 else "other/shift")
print("injective window-8 maps:", names)
# check of the injectivity test: a sample of non-injective conserving maps identify two periodic points, the injective ones do not
import itertools
def apply_cyc(f, w, bits):
    n = len(bits); out = []
    for i in range(n):
        e = 0
        for k in range(w): e = (e << 1) | bits[(i + k) % n]
        out.append(f[e])
    return tuple(out)
bad = 0; noninj_witness = 0; tested = 0; nnon = 0
import random; random.seed(3)
for f in random.sample(sols, 150) + inj:
    isinj = injective(8, f)
    seen = {}
    coll = False
    for n in range(1, 14):
        for m in range(1 << n):
            if (m & (m >> 1)) or ((m & 1) and (m >> (n - 1)) & 1): continue
            bits = tuple((m >> (n - 1 - i)) & 1 for i in range(n))
            im = apply_cyc(f, 8, bits)
            if (n, im) in seen and seen[(n, im)] != bits: coll = True
            seen[(n, im)] = bits
    tested += 1; nnon += (not isinj)
    if isinj and coll: bad += 1
    if (not isinj) and coll: noninj_witness += 1
print("sample of %d window-8 maps: injective-by-graph but colliding on periodic points: %d; non-injective with an explicit collision of period <=13: %d of %d non-injective" % (tested, bad, noninj_witness, nnon), flush=True)
