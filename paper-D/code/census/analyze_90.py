"""90R, radius R: structure of the full solution space (from linear_90.kernel_basis).
 (a) dimension of the affine part = dim( solutions  ∩  {tables of affine maps} ), computed algebraically;
 (b) algebraic normal form (ANF) of the output bits of a basis of a complement of the affine part;
 (c) enumeration of ALL 2^dim solutions with cheap tags: named sigma^a F^b, constant, Hedlund-balanced
     (each output symbol exactly 4^{2R} times -- necessary for surjectivity, hence for invertibility);
     pair-graph injectivity test only on the balanced ones. Layer structure on all of them.
"""
import sys, time, json, itertools
from collections import Counter
from centralizer import *
from linear_90 import build_rows, kernel_basis, vec_to_table

def affine_table_basis(R):
    """tables (as bit vectors, same layout as linear_90) of the affine maps: constants and the 2(2R+1) input bits,
    for each of the 2 output bits."""
    n = 2 * R + 1; N = 4 ** n; vecs = []
    for ob in (0, 1):                       # output bit: 0 -> X (vector bit 2e), 1 -> Y (bit 2e+1)
        pos = lambda e: 2 * e + (0 if ob == 0 else 1)
        vecs.append(sum(1 << pos(e) for e in range(N)))                 # constant 1
        for cell in range(n):
            for ib in (1, 0):                                             # input x-bit (1) or y-bit (0) of that cell
                vecs.append(sum(1 << pos(e) for e in range(N) if (index_word(e, n)[cell] >> ib) & 1))
    return vecs

def rank_of(vecs):
    piv = {}
    for r in vecs:
        while r:
            p = r & -r
            if p in piv: r ^= piv[p]
            else: piv[p] = r; break
    return len(piv)

def anf(table_bit, n_in):
    """Moebius transform: coefficients of the ANF of a boolean function given as a list indexed by input int."""
    a = list(table_bit); N = len(a)
    step = 1
    while step < N:
        for i in range(N):
            if i & step: a[i] ^= a[i ^ step]
        step <<= 1
    return [i for i in range(N) if a[i]]

def monomial_name(mask, R):
    n = 2 * R + 1; names = []
    for cell in range(n):
        for ib, nm in ((1, "x"), (0, "y")):
            bit = 2 * (n - 1 - cell) + ib
            if (mask >> bit) & 1: names.append(f"{nm}{cell - R:+d}")
    return "*".join(names) if names else "1"

def run(R, enumerate_all):
    t0 = time.time()
    rows, nvar = build_rows(R); rank, basis = kernel_basis(rows, nvar); d = len(basis)
    aff = affine_table_basis(R)
    ra = rank_of(aff); ru = rank_of(aff + basis)
    d_aff = ra + d - ru
    print(f"=== 90R R={R}: solution space dim {d} (2^{d} = {2**d} maps); affine-map table space dim {ra}; "
          f"dim of affine solutions = {d_aff} (2^{d_aff} = {2**d_aff}); non-affine complement dim {d - d_aff}  [{time.time()-t0:.1f}s]")
    # complement basis: reduce solution basis against the affine space
    piv = {}
    for r in aff:
        while r:
            p = r & -r
            if p in piv: r ^= piv[p]
            else: piv[p] = r; break
    comp = []
    for v in basis:
        r = v
        while r:
            p = r & -r
            if p in piv: r ^= piv[p]
            else: piv[p] = r; comp.append(v); break
    print(f"  {len(comp)} basis vectors lie outside the affine space; ANF of their output bits (inputs x_j,y_j, j = offset):")
    n = 2 * R + 1; f, F_loc, Finv_loc = make_rule(90)
    for v in comp[:6]:
        tab = vec_to_table(v, R)
        g = BlockMap(R, tab)
        xs = anf([(tab[e] >> 1) for e in range(4 ** n)], 2 * n); ys = anf([(tab[e] & 1) for e in range(4 ** n)], 2 * n)
        deg = max([bin(m).count("1") for m in xs + ys] + [0])
        print(f"   X_out = " + " + ".join(monomial_name(m, R) for m in xs)[:300])
        print(f"   Y_out = " + " + ".join(monomial_name(m, R) for m in ys)[:300] + f"      (degree {deg}, balanced={Counter(tab).most_common(1)[0][1] == 4**(2*R)})")
    if not enumerate_all: return
    names = named_maps(R, F_loc, Finv_loc)
    counts = Counter(); N = 4 ** n; inv_list = []
    t0 = time.time()
    tables_by_bal = []
    for mask in range(1 << d):
        v = 0
        for i in range(d):
            if mask >> i & 1: v ^= basis[i]
        tab = vec_to_table(v, R)
        key = tuple(tab)
        c = Counter(tab); balanced = len(c) == 4 and all(x == 4 ** (2 * R) for x in c.values())
        name = names.get(key)
        if name is not None: counts["named"] += 1
        elif len(c) == 1: counts["const"] += 1
        else: counts["other"] += 1
        counts["balanced"] += balanced
        if balanced:
            g = BlockMap(R, tab)
            inj = True if name is not None else injective_pairgraph_np(g)
            if inj:
                counts["invertible"] += 1
                counts["invertible_affine"] += is_affine(g)
                inv_list.append({"label": f"sigma^{name[0]} F^{name[1]}" if name else "OTHER", "affine": is_affine(g),
                                 "layer": layer_structure(g)[0], "table": g.table_str()})
        if mask % 8192 == 0: print(f"    ... {mask}/{1<<d} [{time.time()-t0:.0f}s] {dict(counts)}", flush=True)
    print(f"  enumeration of all 2^{d} maps done [{time.time()-t0:.0f}s]: {dict(counts)}")
    print(f"  invertible commuting maps: {counts['invertible']} (affine: {counts['invertible_affine']}, non-affine: {counts['invertible']-counts['invertible_affine']})")
    for r in inv_list:
        print(f"    {r['label']:14s} affine={r['affine']} {r['layer']:15s} table={r['table'] if R <= 1 else r['table'][:64]+'...'}")
    json.dump({"R": R, "dim": d, "dim_affine": d_aff, "counts": dict(counts), "invertible": inv_list},
              open(f"results_analyze90_R{R}.json", "w"), indent=1)

if __name__ == "__main__":
    R = int(sys.argv[1]); run(R, enumerate_all=(len(sys.argv) > 2 and sys.argv[2] == "all"))
