"""Radius 2: compares the two implementations of (T_P) (scan_r2.TPtab, checks3.TP) on 153
functions g for P = 1, 2, and counts the non-affine g, those with (N), those with (N) and (D),
and the self-dual non-affine g."""
import random
from scan_r2 import TPtab, is_affine
from checks3 import TP, J_of, selfdual
random.seed(1)
bad = 0
for tab in random.sample(range(65536), 150) + [0xff00, 0x00ff, 0x0339]:
    g = lambda x, tab=tab: (tab >> (x[0] | (x[1] << 1) | (x[2] << 2) | (x[3] << 3))) & 1
    for P in (1, 2):
        if TP(g, 2, P)[0] != TPtab(tab, P): bad += 1
print("disagreements between the two T_P implementations:", bad)
cnt = {"nonaffine": 0, "N": 0, "N&D": 0, "sd_nonaffine": 0}
for tab in range(65536):
    if is_affine(tab): continue
    g = lambda x, tab=tab: (tab >> (x[0] | (x[1] << 1) | (x[2] << 2) | (x[3] << 3))) & 1
    cnt["nonaffine"] += 1
    J, const = J_of(g, 2)
    sd = selfdual(g, 2)
    cnt["sd_nonaffine"] += sd
    if (not const) and J is not None and J >= 2:
        cnt["N"] += 1
        if not sd: cnt["N&D"] += 1
print(cnt)
