"""Psi_m(sigma): in every 'train' -- a maximal configuration piece  W s_1 33 s_2 33 ... 33 s_m W  where each s_j != 3,
the inner walls are exactly '33', and the outer walls W are runs of >= 3 threes -- replace (s_1..s_m) by
(s_{sigma^-1(1)}, ..., s_{sigma^-1(m)}). Check on random cyclic configurations: commutes with F, is a bijection
(sigma^-1 inverts), and sigma -> Psi_m(sigma) is a homomorphism."""
import random, itertools
from centralizer import make_rule
f, F_loc, _ = make_rule(54)
def F(c):
    L = len(c); return [F_loc(c[i-1], c[i], c[(i+1) % L]) for i in range(L)]
def trains(c):
    """list of (positions of the m cars) for each train in cyclic c."""
    L = len(c); out = []
    for i in range(L):
        # train starts at car position i: preceded by >= 3 threes
        if c[i] != 3 and all(c[(i-k) % L] == 3 for k in (1, 2, 3)):
            pos = [i]; j = i
            while True:
                if c[(j+1) % L] == 3 and c[(j+2) % L] == 3 and c[(j+3) % L] != 3 and c[(j+4) % L] == 3 and c[(j+5) % L] == 3 and not (c[(j+4) % L] == 3 and c[(j+5) % L] == 3 and c[(j+6) % L] == 3 and False):
                    # next car at j+3, inner wall exactly two threes (j+1, j+2) and car followed by a wall
                    pass
                break
            # rebuild properly: walk cars
            pos = [i]; j = i
            while True:
                if not (c[(j+1) % L] == 3 and c[(j+2) % L] == 3): pos = None; break      # car must be followed by a wall
                if c[(j+3) % L] == 3:                                                  # wall of length >= 3: train ends
                    break
                pos.append((j+3) % L); j = (j+3) % L
                if len(pos) > L: pos = None; break
            if pos: out.append(pos)
    return out
def Psi(c, m, sigma):
    d = list(c)
    for pos in trains(c):
        if len(pos) == m:
            vals = [c[p] for p in pos]
            for k in range(m): d[pos[sigma[k]]] = vals[k]
    return d
rng = random.Random(0)
def rand_cfg():
    # concatenate random blocks: walls of length 2-4, cars, longer segments
    c = []
    while len(c) < 60:
        r = rng.random()
        c += [3] * rng.choice([2, 2, 3, 4]); c += [rng.choice([0, 1, 2])] if r < 0.7 else [rng.choice([0,1,2]) for _ in range(rng.randint(2, 4))]
    return c
ok_comm = ok_bij = ok_hom = True; nontriv = 0
for trial in range(3000):
    m = rng.randint(2, 4); c = rand_cfg()
    s1 = list(range(m)); rng.shuffle(s1); s2 = list(range(m)); rng.shuffle(s2)
    a = Psi(F(c), m, s1); b = F(Psi(c, m, s1))
    if a != b: ok_comm = False; print("commutation fails", m, c); break
    inv = [0] * m
    for k in range(m): inv[s1[k]] = k
    if Psi(Psi(c, m, s1), m, inv) != c: ok_bij = False; break
    comp = [s1[s2[k]] for k in range(m)]
    if Psi(Psi(c, m, s2), m, s1) != Psi(c, m, comp): ok_hom = False; break
    nontriv += Psi(c, m, s1) != c
print(f"commutes with F: {ok_comm}; sigma^-1 inverts: {ok_bij}; homomorphism: {ok_hom}; nontrivial on {nontriv}/3000 samples")
