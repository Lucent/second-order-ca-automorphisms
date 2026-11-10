"""The non-affine radius-2 rules with no (T_P), P<=3, read from out_scan_r2.txt:
test (T_4), and (T_5) for those failing (T_4), within a time limit of 480 s."""
import re, sys, time
from scan_r2 import TPtab
tabs = [int(m.group(1)) for m in re.finditer(r"tab=(\d+) .*affine=False", open("out_scan_r2.txt").read())]
t0 = time.time(); res = {}
for P in (4, 5):
    for tab in tabs:
        if res.get(tab): continue
        if time.time() - t0 > 480: break
        if TPtab(tab, P): res[tab] = P
    print("after P=%d: certified %d of %d; remaining: %s" % (P, len(res), len(tabs), ["0x%04x" % t for t in tabs if t not in res]), flush=True)
print({("0x%04x" % t): p for t, p in res.items()})
