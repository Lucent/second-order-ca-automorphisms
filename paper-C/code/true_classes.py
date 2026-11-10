"""True symmetry classes of second-order ECA: reflection and kappa: f -> f(1-a,1-b,1-c) (no outer complement).
Check conjugacy claim numerically, recount walled rules per class."""
import itertools, json
from centralizer import make_rule
bit = lambda r, a, b, c: (r >> ((a << 2) | (b << 1) | c)) & 1
mk = lambda g: sum(g(a, b, c) << ((a << 2) | (b << 1) | c) for a in (0,1) for b in (0,1) for c in (0,1))
mir = lambda r: mk(lambda a, b, c: bit(r, c, b, a)); kap = lambda r: mk(lambda a, b, c: bit(r, 1-a, 1-b, 1-c))
# conjugacy check: kappa o F_f = F_{kap f} o kappa on all 3-windows
for r in (110, 30, 54, 37):
    _, F1, _ = make_rule(r); _, F2, _ = make_rule(kap(r))
    ok = all(3 - F1(a, b, c) == F2(3 - a, 3 - b, 3 - c) for a in range(4) for b in range(4) for c in range(4))
    print(f"kappa conjugates {r}R to {kap(r)}R: {ok}")
def wall(r):
    _, F, _ = make_rule(r)
    for k in (1, 2, 3):
        for w in itertools.product(range(4), repeat=k):
            if all(tuple(F(*ctx[i:i+3]) for i in range(k)) == w for ctx in (tuple([a]) + w + tuple([b]) for a in range(4) for b in range(4))):
                return "".join(map(str, w))
seen = set(); classes = []
for r in range(256):
    if r in seen: continue
    cl = sorted({r, mir(r), kap(r), mir(kap(r))}); seen |= set(cl); classes.append(cl)
walls = {r: wall(r) for r in range(256)}
const = all(len({bool(walls[x]) for x in cl}) == 1 for cl in classes)
wc = [cl for cl in classes if walls[cl[0]]]; nw = [cl for cl in classes if not walls[cl[0]]]
print(f"{len(classes)} true classes; walled rules {sum(1 for r in walls if walls[r])}; walled classes {len(wc)}; wall-free classes {len(nw)}; property constant on classes: {const}")
print("walled reps:", " ".join(f"{cl[0]}:{walls[cl[0]]}" for cl in wc))
print("wall-free reps:", " ".join(str(cl[0]) for cl in nw))
print("class of 110:", [cl for cl in classes if 110 in cl], " class of 30:", [cl for cl in classes if 30 in cl], " class of 54:", [cl for cl in classes if 54 in cl], " class of 37:", [cl for cl in classes if 37 in cl])
json.dump({"walled": wc, "wallfree": nw, "walls": walls}, open("true_classes.json", "w"))
