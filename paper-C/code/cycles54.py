import itertools
from collections import Counter
from centralizer import make_rule
f, F_loc, _ = make_rule(54)
def F_seg(s):
    w = (3,) + s + (3,); return tuple(F_loc(*w[i:i+3]) for i in range(len(s)))
for L in range(1, 11):
    S = [s for s in itertools.product(range(4), repeat=L) if s[0] != 3 and s[-1] != 3 and not any(s[i] == 3 and s[i+1] == 3 for i in range(L-1))]
    seen = set(); lens = Counter()
    for s in S:
        if s in seen: continue
        k = 0; t = s
        while True:
            seen.add(t); t = F_seg(t); k += 1
            if t == s: break
        lens[k] += 1
    kbest, mbest = max(lens.items(), key=lambda kv: kv[1])
    print(f"L={L}: |S_L|={len(S)}, cycles={sum(lens.values())}, distinct lengths={len(lens)}, max length={max(lens)}, most repeated: {mbest} cycles of length {kbest}", flush=True)
