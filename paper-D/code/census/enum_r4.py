"""110R radius 4: enumerate with minisat22, dump tables immediately, classify with numpy checks (fast)."""
import sys, time, json
from centralizer import *
from classify_r4 import check_commutes_np_R
rule, R = int(sys.argv[1]), 4
f, F_loc, Finv_loc = make_rule(rule)
sols, capped, nv, nc, t_enc, t_enum = enumerate_solutions(R, f, F_loc, 1000, "minisat22")
print(f"=== rule {rule}R radius {R} [minisat22]: {nv} variables, {nc} clauses, encode {t_enc:.1f}s, enumeration {t_enum:.1f}s, {len(sols)} solutions", flush=True)
with open(f"raw_solutions_{rule}_R4.txt", "w") as fh:
    for t in sols: fh.write("".join(map(str, t)) + "\n")
t0 = time.time(); names = named_maps(R, F_loc, Finv_loc); print(f"  {len(names)} named maps built [{time.time()-t0:.0f}s]", flush=True)
counts = {"named": 0, "const": 0, "other": 0, "layerwise-same": 0, "joint": 0, "layer-separable": 0, "cross": 0, "commute_ok": 0}
rows = []
for tab in sols:
    g = BlockMap(R, tab); counts["commute_ok"] += check_commutes_np_R(tab, F_loc, R)
    name = names.get(g.key())
    if name: label = f"sigma^{name[0]} F^{name[1]}"; counts["named"] += 1
    elif len(set(tab)) == 1: label = f"const ({tab[0]>>1},{tab[0]&1})"; counts["const"] += 1
    else: label = "OTHER"; counts["other"] += 1
    k = layer_structure(g)[0]; counts[k] += 1; rows.append((label, k))
missing = [v for kk, v in names.items() if kk not in {tuple(t) for t in sols}]
print(f"  classification [{time.time()-t0:.0f}s]: {counts}; expected named {len(names)}; all named found: {not missing}", flush=True)
for label, k in sorted(rows): print(f"    {label:16s} {k}")
json.dump({"rule": rule, "R": R, "nvars": nv, "nclauses": nc, "t_encode": t_enc, "t_enum": t_enum, "nsol": len(sols),
           "counts": counts, "solutions": rows}, open(f"results_R4_{rule}_fast.json", "w"), indent=1)
