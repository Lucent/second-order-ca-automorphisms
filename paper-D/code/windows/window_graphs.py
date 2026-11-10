"""Window graphs W^{(k,p)}_M of the second-order automaton F_f of an elementary rule f.
Vertex: word omega in A^{M+1} on [q, q+M].  Edge: for every choice of the outside cells needed,
new vertex = (sigma^k F^p)(...) restricted to [q, q+M], i.e. F^p-image on [q+k, q+k+M].
(k,p)=(0,1): zero-speed window graph V_M (two-sided input).  (k,p)=(1,1): right light-speed graph
(input only on the right), (-1,1): left light-speed graph.
graph(rule, M, k, p) returns (number of SCCs, number of terminal SCCs, sizes of the largest terminal SCCs);
the graph is strongly connected if and only if the first number is 1.

Run as a script it prints, for one representative of each of the 88 classes of rules under reflection
and the cellwise complement, one line with three strings: for (k,p) = (0,1), (1,1), (-1,1) and each M
in the list (default 0,...,6), the letter S if W^{(k,p)}_M is strongly connected and a dot if not.

Usage: python window_graphs.py [M1,M2,...] [first_rule last_rule]"""
import numpy as np, json, sys

def scc(nv, src, tgt):
    """iterative Tarjan; returns (ncomp, labels)"""
    import sys
    order = np.argsort(src, kind='stable'); s2 = src[order]; t2 = tgt[order]
    ptr = np.searchsorted(s2, np.arange(nv+1))
    adj = [np.unique(t2[ptr[v]:ptr[v+1]]).tolist() for v in range(nv)]
    index = [-1]*nv; low = [0]*nv; onst = [False]*nv; st = []; lab = [-1]*nv
    idx = 0; nc = 0
    for root in range(nv):
        if index[root] != -1: continue
        stack = [(root, 0)]
        index[root] = low[root] = idx; idx += 1; st.append(root); onst[root] = True
        while stack:
            v, i = stack[-1]
            if i < len(adj[v]):
                stack[-1] = (v, i+1); w = adj[v][i]
                if index[w] == -1:
                    index[w] = low[w] = idx; idx += 1; st.append(w); onst[w] = True
                    stack.append((w, 0))
                elif onst[w]:
                    low[v] = min(low[v], index[w])
            else:
                stack.pop()
                if stack: low[stack[-1][0]] = min(low[stack[-1][0]], low[v])
                if low[v] == index[v]:
                    while True:
                        w = st.pop(); onst[w] = False; lab[w] = nc
                        if w == v: break
                    nc += 1
    return nc, np.array(lab)

def Fstep(rule, W):
    """W: (n_words, L) int array of cells 0..3 -> (n_words, L-2)"""
    x = W >> 1; y = W & 1
    a, b, c = y[:, :-2], y[:, 1:-1], y[:, 2:]
    fv = (rule >> (4*a + 2*b + c)) & 1
    return (b << 1) | (x[:, 1:-1] ^ fv)

def graph(rule, M, k, p, inv=False):
    n = M + 1
    lo = min(0, k - p); hi = max(0, k + p)        # cells needed: [q+k-p, q+k+M+p]
    left = max(0, p - k)                           # number of cells needed left of q
    right = max(0, k + p)                          # number needed right of q+M
    # cells needed: q+k-p .. q+k+M+p ; outside window: those < q or > q+M
    L = n + left + right
    nv = 4 ** n; ni = 4 ** (left + right)
    V = np.arange(nv)
    vd = np.stack([(V >> (2*(n-1-j))) & 3 for j in range(n)], axis=1)
    I = np.arange(ni)
    idg = np.stack([(I >> (2*(left+right-1-j))) & 3 for j in range(left+right)], axis=1) if left+right else np.zeros((1,0),int)
    W = np.empty((nv * ni, L), dtype=np.int64)
    vv = np.repeat(vd, ni, axis=0); ii = np.tile(idg, (nv, 1))
    W[:, :left] = ii[:, :left]; W[:, left:left+n] = vv; W[:, left+n:] = ii[:, left:]
    # W covers [q-left, q+M+right]; we need F^p on [q+k-p, q+k+M+p]
    start = (k - p) + left
    W = W[:, start:start + n + 2*p]
    rr = rule
    for _ in range(p):
        W = Fstep(rr, W) if not inv else Finvstep(rr, W)
    tgt = np.zeros(len(W), dtype=np.int64)
    for j in range(n): tgt = (tgt << 2) | W[:, j]
    src = np.repeat(V, ni)
    ncomp, lab = scc(nv, src, tgt)
    # terminal SCCs: no edge leaving
    out_other = np.zeros(ncomp, bool)
    m = lab[src] != lab[tgt]
    out_other[np.unique(lab[src][m])] = True
    term = [c for c in range(ncomp) if not out_other[c]]
    sizes = np.bincount(lab, minlength=ncomp)
    return ncomp, len(term), sorted([int(sizes[c]) for c in term], reverse=True)[:6]

def Finvstep(rule, W):
    x = W >> 1; y = W & 1
    a, b, c = x[:, :-2], x[:, 1:-1], x[:, 2:]
    fv = (rule >> (4*a + 2*b + c)) & 1
    return ((y[:, 1:-1] ^ fv) << 1) | b

if __name__ == '__main__':
    def bit(r, a, b, c): return (r >> ((a << 2) | (b << 1) | c)) & 1
    def mk(g): return sum(g(a, b, c) << ((a << 2) | (b << 1) | c) for a in (0, 1) for b in (0, 1) for c in (0, 1))
    mirror = lambda r: mk(lambda a, b, c: bit(r, c, b, a)); kappa = lambda r: mk(lambda a, b, c: bit(r, 1 - a, 1 - b, 1 - c))
    reps = sorted({min(r, mirror(r), kappa(r), mirror(kappa(r))) for r in range(256)})
    Ms = [int(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0, 1, 2, 3, 4, 5, 6]
    lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, 255)
    for rule in reps:
        if not lo <= rule <= hi: continue
        sc = lambda k, p: ''.join('S' if graph(rule, M, k, p)[0] == 1 else '.' for M in Ms)
        print(f"{rule:4d}  zero-speed (0,1): {sc(0, 1)}   (1,1): {sc(1, 1)}   (-1,1): {sc(-1, 1)}", flush=True)
