"""Exact centralizer of 90R (f = y_{i-1} xor y_{i+1}) among ALL block maps of radius R, by GF(2) Gaussian
elimination. Because f is linear, the commutation constraints of centralizer.py
    X_eL = Y_e2 ,   Y_eL xor X_e2 = Y_e1 xor Y_e3
are homogeneous F_2-linear equations in the 2*4^(2R+1) table bits, so the set of commuting block maps
(of ANY form, not only linear ones) is an F_2-subspace of table space; its size is 2^dim(kernel).
Theory: the linear commutant of M = [[0,1],[1,p]] (p = s^-1 + s) in M_2(F_2[s,s^-1]) is {aI + bM};
with radius R, a has 2R+1 free coefficients and b has 2R-1, and the constant map to the fixed point
(1,1) adds one more dimension: predicted dim = 4R+1.  Invertible affine maps: det(aI+bM) = a^2+abp+b^2
must be a unit (a monomial).
"""
import sys, time, json, itertools, random
from centralizer import *

def build_rows(R):
    n = 2 * R + 1; m = 2 * R + 3
    f, F_loc, _ = make_rule(90)
    X = lambda e: 1 << (2 * e); Y = lambda e: 1 << (2 * e + 1)
    rows = set()
    for ew in range(4 ** m):
        w = index_word(ew, m)
        eL = word_index(apply_local(F_loc, w, 1))
        e1, e2, e3 = word_index(w[0:n]), word_index(w[1:n + 1]), word_index(w[2:n + 2])
        r1 = X(eL) ^ Y(e2); r2 = Y(eL) ^ X(e2) ^ Y(e1) ^ Y(e3)
        if r1: rows.add(r1)
        if r2: rows.add(r2)
    return rows, 2 * 4 ** n

def kernel_basis(rows, nvar):
    piv = {}
    for r in rows:
        while r:
            p = r & -r
            if p in piv: r ^= piv[p]
            else: piv[p] = r; break
    # reduced echelon: eliminate pivot bits from all other rows
    for p in sorted(piv, reverse=True):
        r = piv[p]
        for q in piv:
            if q != p and piv[q] & p: piv[q] ^= r
    rank = len(piv)
    free = [i for i in range(nvar) if (1 << i) not in piv]
    basis = []
    for fv in free:                      # set free var fv = 1, others 0, solve pivots
        vec = 1 << fv
        for p, r in piv.items():
            if r & (1 << fv): vec |= p
        basis.append(vec)
    return rank, basis

def vec_to_table(vec, R):
    n = 2 * R + 1
    return [(((vec >> (2 * e)) & 1) << 1) | ((vec >> (2 * e + 1)) & 1) for e in range(4 ** n)]

# ---- Laurent polynomial arithmetic over F_2 (set of exponents)
def pmul(a, b):
    out = set()
    for i in a:
        for j in b: out ^= {i + j}
    return out
def padd(a, b): return set(a) ^ set(b)
P = {-1, 1}
def count_invertible_affine(R):
    """(a,b) with a supported in [-R,R], b in [-(R-1),R-1], det = a^2 + a b p + b^2 a monomial."""
    cnt = 0; total = 0
    ea = list(range(-R, R + 1)); eb = list(range(-(R - 1), R))
    for ma in range(1 << len(ea)):
        a = {ea[i] for i in range(len(ea)) if ma >> i & 1}
        for mb in range(1 << len(eb)):
            b = {eb[i] for i in range(len(eb)) if mb >> i & 1}
            total += 1
            det = padd(padd(pmul(a, a), pmul(pmul(a, b), P)), pmul(b, b))
            if len(det) == 1: cnt += 1
    return cnt, total

def run(R, classify_all):
    t0 = time.time()
    rows, nvar = build_rows(R)
    t_rows = time.time() - t0; t0 = time.time()
    rank, basis = kernel_basis(rows, nvar)
    t_elim = time.time() - t0
    d = len(basis)
    print(f"\n=== 90R radius {R}: {nvar} unknowns, {len(rows)} distinct equations, rank {rank}, kernel dim {d} "
          f"=> exactly 2^{d} = {2**d} commuting block maps (predicted 4R+1 = {4*R+1}); rows {t_rows:.1f}s, elimination {t_elim:.1f}s", flush=True)
    f, F_loc, Finv_loc = make_rule(90)
    names = named_maps(R, F_loc, Finv_loc)
    # basis vectors: independent commutation check + affinity (sum of affine maps is affine, so this covers everything)
    t0 = time.time()
    bas_aff = 0
    for v in basis:
        g = BlockMap(R, vec_to_table(v, R))
        assert check_commutes(g, F_loc, R), "kernel vector does not commute"
        bas_aff += is_affine(g)
    print(f"  all {d} basis maps commute (independent re-check); {bas_aff}/{d} affine => all 2^{d} solutions are affine over F_2"
          if bas_aff == d else f"  WARNING: only {bas_aff}/{d} basis maps affine")
    ninv_th, tot = count_invertible_affine(R)
    print(f"  theory: linear maps aI+bM of radius {R}: {tot}; invertible among them (det a unit): {ninv_th}; "
          f"with the constant (1,1) added: {2*tot} maps, {2*ninv_th} invertible")
    out = {"R": R, "nvar": nvar, "neq": len(rows), "rank": rank, "dim": d, "count": 2 ** d, "t_rows": t_rows, "t_elim": t_elim,
           "basis_all_affine": bas_aff == d, "theory_linear": tot, "theory_invertible": 2 * ninv_th}
    if classify_all:
        counts = {"named": 0, "const": 0, "other": 0, "invertible": 0, "layerwise-same": 0, "layer-separable": 0, "cross": 0, "joint": 0, "affine": 0}
        sols = []
        for mask in range(1 << d):
            v = 0
            for i in range(d):
                if mask >> i & 1: v ^= basis[i]
            g = BlockMap(R, vec_to_table(v, R))
            assert check_commutes(g, F_loc, R)
            name = names.get(g.key()); tab = g.table
            if name is not None: label = f"sigma^{name[0]} F^{name[1]}"; counts["named"] += 1
            elif len(set(tab)) == 1: label = f"const ({tab[0]>>1},{tab[0]&1})"; counts["const"] += 1
            else: label = "OTHER"; counts["other"] += 1
            inv = (True if name is not None else False if len(set(tab)) == 1 else injective_pairgraph_np(g))
            counts["invertible"] += inv
            kind, _ = layer_structure(g); counts[kind] += 1
            aff = is_affine(g); counts["affine"] += aff
            sols.append({"label": label, "invertible": inv, "layer": kind, "affine": aff, "table": g.table_str() if R <= 2 else ""})
        print(f"  full classification of all {2**d} maps ({time.time()-t0:.1f}s): {counts}; named expected {len(names)}")
        out["counts"] = counts; out["solutions"] = sols
    return out

def rank_only(R):
    """streaming GF(2) rank of the commutation system (no basis), for large R."""
    n = 2 * R + 1; m = 2 * R + 3; nvar = 2 * 4 ** n
    f, F_loc, _ = make_rule(90)
    X = lambda e: 1 << (2 * e); Y = lambda e: 1 << (2 * e + 1)
    piv = {}; t0 = time.time(); nrows = 0
    for ew in range(4 ** m):
        w = index_word(ew, m)
        eL = word_index(apply_local(F_loc, w, 1))
        e1, e2, e3 = word_index(w[0:n]), word_index(w[1:n + 1]), word_index(w[2:n + 2])
        for r in (X(eL) ^ Y(e2), Y(eL) ^ X(e2) ^ Y(e1) ^ Y(e3)):
            nrows += 1
            while r:
                p = r & -r
                if p in piv: r ^= piv[p]
                else: piv[p] = r; break
        if ew % 50000 == 0: print(f"    ... {ew}/{4**m} windows, rank so far {len(piv)} [{time.time()-t0:.0f}s]", flush=True)
    d = nvar - len(piv)
    print(f"=== 90R radius {R}: {nvar} unknowns, {nrows} equations, rank {len(piv)}, kernel dim {d} => exactly 2^{d} commuting maps "
          f"(prediction 2^(2R)+1 = {4**R+1}); {time.time()-t0:.0f}s", flush=True)
    return d

if __name__ == "__main__":
    if sys.argv[1] == "rank":
        rank_only(int(sys.argv[2])); sys.exit()
    Rs = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [1, 2, 3]
    res = []
    for R in Rs:
        res.append(run(R, classify_all=(R <= 2)))
        json.dump(res, open(f"results_linear90.json", "w"), indent=1)
