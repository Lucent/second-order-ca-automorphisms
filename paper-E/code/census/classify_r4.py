"""numpy commutation check g o F = F o g on all 4^{2R+3} windows, any R (independent of the SAT/propagation encodings)."""
import numpy as np

def check_commutes_np_R(tab, F_loc, R):
    tab = np.asarray(tab, dtype=np.int64); n = 2 * R + 1; m = 2 * R + 3
    idx = np.arange(4 ** m, dtype=np.int64)
    cells = [(idx >> (2 * (m - 1 - i))) & 3 for i in range(m)]
    Ft = np.array([[[F_loc(a, b, c) for c in range(4)] for b in range(4)] for a in range(4)], dtype=np.int64)
    Fw = [Ft[cells[i], cells[i + 1], cells[i + 2]] for i in range(m - 2)]
    eL = sum(Fw[i] << (2 * (n - 1 - i)) for i in range(n)); lhs = tab[eL]
    e = [sum(cells[j + i] << (2 * (n - 1 - i)) for i in range(n)) for j in range(3)]
    g = [tab[x] for x in e]; rhs = Ft[g[0], g[1], g[2]]
    return bool((lhs == rhs).all())
