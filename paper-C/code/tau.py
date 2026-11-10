"""The wall-trading involution tau for 54R:
   3 3 3 s 3 3   <->   3 3 pi(s) s 3 3        (s in {0,1,2}, pi = 0->1->2->0)
i.e. a singleton segment whose left wall has length >= 3 becomes the two-cell segment (pi(s), s), and a two-cell
segment of the form (pi(s), s) becomes a singleton with a longer left wall. Everything else is fixed."""
import itertools, random, json
import numpy as np
from centralizer import make_rule, index_word, word_index, BlockMap, injective_pairgraph_np, injective_periodic, named_maps
from classify_r4 import check_commutes_np_R
f, F_loc, Finv_loc = make_rule(54); R = 3; n = 7; N = 4 ** n
PI = {0: 1, 1: 2, 2: 0, 3: 3}
def tau_local(w):                      # w: 7 cells, centre index 3
    c = w[3]
    # centre is the wall cell just left of a singleton with left wall >= 3:  . 3 3 [3] s 3 3
    if c == 3 and w[1] == 3 and w[2] == 3 and w[4] != 3 and w[5] == 3 and w[6] == 3: return PI[w[4]]
    # centre is a = pi(s) in a two-cell segment (a, s):  . 3 3 [a] s 3 3  with a = pi(s)
    if c != 3 and w[1] == 3 and w[2] == 3 and w[4] != 3 and w[5] == 3 and w[6] == 3 and c == PI[w[4]]: return 3
    return c
tab = [tau_local(index_word(e, n)) for e in range(N)]
ident = [index_word(e, n)[3] for e in range(N)]
print("changed entries:", sum(1 for a, b in zip(tab, ident) if a != b))
print("commutes with F:", check_commutes_np_R(tab, F_loc, R))
g = BlockMap(R, tab)
print("injective (pair graph):", injective_pairgraph_np(g), "; periodic<=10:", injective_periodic(g, 10))
# involution: tau(tau(c)) = c on random cyclic configurations
def apply(tab, c):
    L = len(c); return [tab[word_index(tuple(c[(i + d) % L] for d in range(-3, 4)))] for i in range(L)]
rng = random.Random(1); ok = True
for _ in range(2000):
    c = [rng.choice([3, 3, 3, 0, 1, 2]) for _ in range(30)]
    if apply(tab, apply(tab, c)) != c: ok = False; break
print("involution on 2000 random cyclic configurations:", ok)
names = named_maps(R, F_loc, Finv_loc); print("tau is a named map:", tuple(tab) in names)
json.dump("".join(map(str, tab)), open("tau_table.json", "w"))
# does tau commute with D = Phi_1(pi)?  compose as radius-6 tables on random configs
def D_local(w):
    return PI[w[3]] if (w[1] == 3 and w[2] == 3 and w[3] != 3 and w[4] == 3 and w[5] == 3) else w[3]
tabD = [D_local(index_word(e, n)) for e in range(N)]
comm = all(apply(tab, apply(tabD, c)) == apply(tabD, apply(tab, c)) for c in ([rng.choice([3,3,3,0,1,2]) for _ in range(30)] for _ in range(500)))
print("tau commutes with D (random test):", comm)
