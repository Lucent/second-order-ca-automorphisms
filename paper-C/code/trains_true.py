"""General train construction for a persistent word W of F_f: a car is a single cell s with W immediately on both
sides (W s W), where the only occurrences of W in W s W are its two ends; a train is a maximal chain
W s_1 W s_2 W ... W s_m W bounded on both sides by WW (two adjacent copies, no gap). Permuting the cars of
m-car trains is tested for commutation with F on random cyclic configurations."""
import itertools, random
from centralizer import make_rule
def persistent_word(F_loc):
    for k in (1, 2, 3):
        for w in itertools.product(range(4), repeat=k):
            if all(tuple(F_loc(*ctx[i:i+3]) for i in range(k)) == w for ctx in (tuple([a]) + w + tuple([b]) for a in range(4) for b in range(4))):
                return list(w)
def occ(c, W):
    L, k = len(c), len(W); return [i for i in range(L) if all(c[(i + j) % L] == W[j] for j in range(k))]
def trains(c, W):
    L, k = len(c), len(W); O = set(occ(c, W)); out = []
    for i in O:
        # train starts at occurrence i if preceded by an adjacent occurrence (WW) and followed by a car
        if (i - k) % L not in O: continue
        pos = []; j = i
        while (j + k + 1) % L in O and (j + k) % L not in O and not any(((j + k) - t) % L in O for t in range(1, k)):
            pos.append((j + k) % L); j = (j + k + 1) % L
            if len(pos) > L: break
        if pos and (j + k) % L in O and len(pos) <= L: out.append(pos)
    return out
def run(rule, trials=1500):
    f, F_loc, _ = make_rule(rule); W = persistent_word(F_loc); k = len(W)
    lin = lambda u: [i for i in range(len(u) - k + 1) if u[i:i + k] == W]
    cars = [s for s in range(4) if lin(W + [s] + W) == [0, k + 1]]
    def F(c):
        L = len(c); return [F_loc(c[i - 1], c[i], c[(i + 1) % L]) for i in range(L)]
    rng = random.Random(rule); ok = True; nontriv = 0
    if not cars: return W, cars, None, 0
    for _ in range(trials):
        c = []
        while len(c) < 70:
            r = rng.random()
            if r < 0.5: c += W + W + sum(([rng.choice(cars)] + W for _ in range(rng.randint(1, 4))), [])
            else: c += [rng.randrange(4) for _ in range(rng.randint(1, 5))]
        m = rng.randint(1, 4); sig = list(range(m)); rng.shuffle(sig)
        def Psi(c):
            d = list(c)
            for pos in trains(c, W):
                if len(pos) == m:
                    vals = [c[p] for p in pos]
                    for t in range(m): d[pos[sig[t]]] = vals[t]
            return d
        if Psi(F(c)) != F(Psi(c)): ok = False; break
        nontriv += Psi(c) != c
    return W, cars, ok, nontriv
bad = []; good = []
import json
for r in [cl[0] for cl in json.load(open("true_classes.json"))["walled"]]:
    f, F_loc, _ = make_rule(r)
    if not persistent_word(F_loc): continue
    W, cars, ok, nt = run(r)
    (good if ok and nt > 0 else bad).append((r, "".join(map(str, W)), cars, ok, nt))
print(f"construction verified (commutes, acts nontrivially) for {len(good)} rules:", ", ".join(f"{r}" for r, *_ in good))
print("not verified:", bad)
