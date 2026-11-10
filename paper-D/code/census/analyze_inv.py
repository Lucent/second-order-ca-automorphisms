"""For each OTHER invertible commuting map: nearest sigma^a F^b (Hamming distance on the table) and the
windows where it differs. Also checks the group structure: whether each OTHER map equals F^b o D_k with
D_k = 'apply pi^k to an isolated symbol s in the pattern 33 s 33' (pi = F's 3-cycle 0->1->2->0 on such cells)."""
import sys, json
from centralizer import *

def describe(res):
    rule, R = res["rule"], res["R"]; n = 2 * R + 1
    f, F_loc, Finv_loc = make_rule(rule)
    names = named_maps(R, F_loc, Finv_loc)
    named_tabs = {v: k for k, v in names.items()}
    others = [r for r in res["solutions"] if r["label"] == "OTHER" and r["invertible"]]
    print(f"=== rule {rule}R R={R}: {len(others)} OTHER invertible maps")
    # candidate generators: D_k (k=1,2) = identity except on windows 3 3 s 3 3 (s != 3) where output = pi^k(s)
    pi = {0: 1, 1: 2, 2: 0, 3: 3}
    def D(k):
        tab = []
        for e in range(4 ** n):
            w = index_word(e, n); s = w[R]
            if R >= 2 and w[R - 2] == w[R - 1] == w[R + 1] == w[R + 2] == 3 and s != 3:
                for _ in range(k): s = pi[s]
            tab.append(s)
        return BlockMap(R, tab)
    gens = {f"D_{k}": D(k) for k in (1, 2)}
    products = {}
    for (a, b), ntab in named_tabs.items():
        for gname, g in gens.items():
            if abs(b) + abs(a) <= R:
                # (sigma^a F^b) o D_k : radius R + R... compare as tables only if the composition has radius <= R;
                # D_k has radius 2, so this needs |a|+|b| == 0 ... instead compose at radius R+|a|+|b| and reduce
                big = BlockMap(R, ntab).compose(g, 2 * R)
                N = 2 * (2 * R) + 1; tab = {}; ok = True
                for e in range(4 ** N):
                    w = index_word(e, N); key = word_index(w[R:N - R])
                    if tab.setdefault(key, big.table[e]) != big.table[e]: ok = False; break
                if ok: products[tuple(tab[k] for k in range(4 ** n))] = f"sigma^{a} F^{b} o {gname}"
    for r in others:
        tab = tuple(int(c) for c in r["table"])
        best = min(named_tabs.items(), key=lambda kv: sum(x != y for x, y in zip(kv[1], tab)))
        (a, b), ntab = best
        diff = [(index_word(e, n), ntab[e], tab[e]) for e in range(4 ** n) if ntab[e] != tab[e]]
        ident = products.get(tab, "not a product sigma^a F^b o D_k")
        print(f"  OTHER = {ident}; nearest named sigma^{a} F^{b} differs on {len(diff)} windows: " +
              ", ".join(f"{''.join(map(str, w))}:{x}->{y}" for w, x, y in diff[:8]) + (" ..." if len(diff) > 8 else ""))

if __name__ == "__main__":
    for path in sys.argv[1:]:
        for res in json.load(open(path)): describe(res)
