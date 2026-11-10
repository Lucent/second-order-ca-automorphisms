"""Checks for the wall mechanism of 54R used in the paper."""
import itertools, json
from centralizer import make_rule, BlockMap, index_word, word_index, check_commutes, named_maps
from classify_r4 import check_commutes_np_R
f, F_loc, Finv_loc = make_rule(54)
# (1) rule 54 identity f(a,b,c) = b xor (a or c)
assert all(f(a,b,c) == (b ^ (a | c)) for a in (0,1) for b in (0,1) for c in (0,1))
# (2) rules with f(.,1,1) = f(1,1,.) = 0 (so 33 is a wall)
wall_rules = [r for r in range(256) if all(((r >> ((a<<2)|(b<<1)|c)) & 1) == 0 for (a,b,c) in [(0,1,1),(1,1,1),(1,1,0)])]
print("rules with f(.11)=f(11.)=0:", len(wall_rules), wall_rules)
# (3) F and F^-1 preserve the set of 33 blocks exactly (33 at (i,i+1) in F(c) iff in c): check on all 4-windows
def F_word(w): return tuple(F_loc(*w[i:i+3]) for i in range(len(w)-2))
ok = True
for w in itertools.product(range(4), repeat=4):
    out = F_word(w)        # positions 1,2 of w
    if (out == (3,3)) != (w[1:3] == (3,3)): ok = False
print("33 blocks preserved exactly by F (all 4-windows):", ok)
# (4) F acts on a walled single cell as the 3-cycle pi, 3 fixed
pi = {s: F_loc(3, s, 3) for s in range(4)}; print("F on 3 s 3:", pi)
# (5) period p_L of the all-zero content in a walled segment of length L, and of all contents (lcm)
from math import lcm
def F_seg(seg):  # one step of F on a walled segment (walls = 3 on both sides)
    w = (3,) + tuple(seg) + (3,)
    return tuple(F_loc(*w[i:i+3]) for i in range(len(seg)))
def period(seg):
    s0 = tuple(seg); s = F_seg(s0); k = 1
    while s != s0: s = F_seg(s); k += 1
    return k
rows = []
for L in range(1, 11):
    pz = period((0,)*L)
    allp = 1; mx = 0
    if L <= 7:
        for seg in itertools.product(range(4), repeat=L):
            p = period(seg); allp = lcm(allp, p); mx = max(mx, p)
    rows.append((L, pz, allp if L <= 7 else None, mx if L <= 7 else None))
    print(f"L={L}: period of 0^L = {pz}; order of F on all 4^L contents = {allp if L<=7 else '?'} (max orbit {mx if L<=7 else '?'})")
# (6) D_k at radius 2 commutes with F; differs from identity on exactly 3 entries; order 3
R = 2; n = 5; N = 4 ** n
def Dk_table(k):
    tab = []
    for e in range(N):
        w = index_word(e, n); s = w[2]
        if w[0] == w[1] == w[3] == w[4] == 3 and s != 3:
            for _ in range(k): s = pi[s]
        tab.append(s)
    return tab
names = named_maps(R, F_loc, Finv_loc)
for k in (1, 2):
    g = BlockMap(R, Dk_table(k))
    ident = BlockMap.shift(R, 0)
    print(f"D_{k}: commutes={check_commutes(g, F_loc, R)}, entries differing from identity={sum(a!=b for a,b in zip(g.table, ident.table))}, in named={g.key() in names}")
