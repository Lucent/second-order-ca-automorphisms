"""Radius 2: rule-by-rule comparison of the two computations of the least P <= 4 with (T_P).

For every g (65536 truth tables) the least P is recomputed with scan_r2.TPtab and written to
out_leastP_r2_scan.txt; this table is compared line by line with out_leastP_r2.txt, the
output of tp_direct.py.  Run tp_direct.py first.
"""
from scan_r2 import TPtab

lines = ["# truth table of g (bit index xi0+2xi1+4xi2+8xi3), least P<=4 with (T_P) ('-' if none)"]
for tab in range(65536):
    lvl = "-"
    for P in (1, 2, 3, 4):
        if TPtab(tab, P):
            lvl = str(P); break
    lines.append("0x%04x %s" % (tab, lvl))
open("out_leastP_r2_scan.txt", "w").write("\n".join(lines) + "\n")
other = open("out_leastP_r2.txt").read().split("\n")
if other and other[-1] == "": other.pop()
diff = [(a, b) for a, b in zip(lines, other) if a != b]
print("lines: %d and %d; differing lines: %d" % (len(lines), len(other), len(diff) + abs(len(lines) - len(other))))
for a, b in diff[:20]:
    print("  ", a, "|", b)
hist = {}
for l in lines[1:]:
    hist[l.split()[1]] = hist.get(l.split()[1], 0) + 1
print("histogram of least P:", dict(sorted(hist.items())))
