"""The 88 classes of elementary rules under reflection and the cellwise complement, and, for each
class: whether it is affine; whether it is permutive in an outer variable; the first word of
length at most 3 that F_f fixes in every context (F(a W b) = W on the cells of W for all cells
a, b), if there is one; and the nontrivial permutations of the four cells that commute with F_f
when applied in every cell (a search over the 23 permutations other than the identity, exact:
commutation is checked on all 4^3 windows of three cells).
Expected totals: 88 classes, 256 rules; 31 classes and 81 rules with a fixed word of length at
most 3; 13 classes with a nontrivial commuting cell permutation.

Usage: python classes88.py
"""
import itertools

def bit(r, a, b, c): return (r >> ((a << 2) | (b << 1) | c)) & 1
def mk(g): return sum(g(a, b, c) << ((a << 2) | (b << 1) | c) for a in (0, 1) for b in (0, 1) for c in (0, 1))
def mirror(r): return mk(lambda a, b, c: bit(r, c, b, a))
def kappa(r): return mk(lambda a, b, c: bit(r, 1 - a, 1 - b, 1 - c))
def F(r, s0, s1, s2):
    x, y = s1 >> 1, s1 & 1
    return (y << 1) | (x ^ bit(r, s0 & 1, y, s2 & 1))
def affine(r):
    z = bit(r, 0, 0, 0)
    return all(bit(r, a, b, c) == z ^ (a & (bit(r, 1, 0, 0) ^ z)) ^ (b & (bit(r, 0, 1, 0) ^ z)) ^ (c & (bit(r, 0, 0, 1) ^ z))
               for a in (0, 1) for b in (0, 1) for c in (0, 1))
def outer_permutive(r):
    return all(bit(r, 0, b, c) != bit(r, 1, b, c) for b in (0, 1) for c in (0, 1)) or \
           all(bit(r, a, b, 0) != bit(r, a, b, 1) for a in (0, 1) for b in (0, 1))
def fixed_word(r):
    for k in (1, 2, 3):
        for w in itertools.product(range(4), repeat=k):
            if all(tuple(F(r, *((a,) + w + (b,))[i:i + 3]) for i in range(k)) == w for a in range(4) for b in range(4)):
                return ''.join(map(str, w))
    return None
def cell_perms(r):
    return [p for p in itertools.permutations(range(4)) if p != (0, 1, 2, 3) and
            all(p[F(r, a, b, c)] == F(r, p[a], p[b], p[c]) for a in range(4) for b in range(4) for c in range(4))]

classes = sorted({tuple(sorted({r, mirror(r), kappa(r), mirror(kappa(r))})) for r in range(256)})
print(len(classes), "classes,", sum(len(c) for c in classes), "rules")
nsym = nw = nwr = 0
for cl in classes:
    r = cl[0]; ps = cell_perms(r); nsym += bool(ps); w = fixed_word(r)
    assert all((fixed_word(q) is None) == (w is None) for q in cl)
    nw += w is not None; nwr += len(cl) * (w is not None)
    tags = (["affine"] if affine(r) else []) + (["outer-permutive"] if outer_permutive(r) else [])
    print(f"{r:3d} {list(cl)} {' '.join(tags)} fixed word: {w if w else 'none'}; cell permutations: {ps if ps else 'none'}")
print(nw, "classes and", nwr, "rules with a fixed word of length at most 3")
print(nsym, "classes with a nontrivial commuting cell permutation")
