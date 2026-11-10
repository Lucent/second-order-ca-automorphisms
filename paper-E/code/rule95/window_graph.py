"""Window graphs of 95R: W^{(k,1)}_M, k = +1, -1 (light speed) and 0 (zero speed), on A^{M+1}: an edge from each word on [0,M] to each of
   its images under sigma^k F, the cells outside ranging over all values.  Strong connectivity for M <= MMAX
   (forward and backward reachability from one vertex).  Control: rule 90, M <= 4 (not strongly connected for k = +1, -1).
   usage: window_graph.py [MMAX]"""
import numpy as np, sys
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
def edges(rule, M, k):
    n = M + 1; N = 4 ** n; idx = np.arange(N * 16, dtype=np.int64)
    w = idx >> 4; ext = idx & 15
    # cells: k=+1: positions 0..M from w, M+1, M+2 from ext; image = F(c)_{1..M+1}.  k=-1: positions -2,-1 from ext; image = F(c)_{-1..M-1}
    cell = {}
    for i in range(n): cell[i] = (w >> (2 * (n - 1 - i))) & 3
    if k == 1: cell[M + 1] = ext >> 2; cell[M + 2] = ext & 3; rng = range(1, M + 2)
    elif k == -1: cell[-2] = ext >> 2; cell[-1] = ext & 3; rng = range(-1, M)
    else: cell[-1] = ext >> 2; cell[M + 1] = ext & 3; rng = range(0, M + 1)
    out = np.zeros_like(idx)
    for i in rng:
        x, y = cell[i] >> 1, cell[i] & 1; a, c = cell[i - 1] & 1, cell[i + 1] & 1
        fy = (rule >> (4 * a + 2 * y + c)) & 1
        out = (out << 2) | ((y << 1) | (x ^ fy))
    return w, out, N
def strongly_connected(src, dst, N):
    for s, d in ((src, dst), (dst, src)):
        order = np.argsort(s, kind="stable"); ss = s[order]; dd = d[order]
        start = np.searchsorted(ss, np.arange(N + 1))
        seen = np.zeros(N, dtype=bool); seen[0] = True; front = np.array([0])
        while len(front):
            nb = np.concatenate([dd[start[v]:start[v + 1]] for v in front]) if len(front) < 64 else dd[np.isin(ss, front)]
            nb = np.unique(nb); nb = nb[~seen[nb]]; seen[nb] = True; front = nb
        if not seen.all(): return False
    return True
for rule in (95, 90):
    for k in (1, -1, 0):
        res = []
        for M in range(0, (MMAX if rule == 95 else 4) + 1):
            s, d, N = edges(rule, M, k); res.append(strongly_connected(s, d, N))
        print("rule %dR, W^(%+d,1)_M strongly connected for M=0..%d:" % (rule, k, len(res) - 1), res, flush=True)
