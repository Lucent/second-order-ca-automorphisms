import itertools
from centralizer import make_rule, index_word
from classify_r4 import check_commutes_np_R
f, F_loc, _ = make_rule(54)
def F_seg(seg):
    w = (3,) + tuple(seg) + (3,); return tuple(F_loc(*w[i:i+3]) for i in range(len(seg)))
def GL_table(L):
    Rr = L + 1; nn = 2 * Rr + 1; tab = []
    for e in range(4 ** nn):
        w = index_word(e, nn); c = Rr; out = w[c]
        for i in range(c - L + 1, c + 1):
            j = i + L - 1
            if i - 2 < 0 or j + 2 >= nn: continue
            seg = w[i:j+1]
            if w[i-2] == w[i-1] == 3 and w[j+1] == w[j+2] == 3 and seg[0] != 3 and seg[-1] != 3 \
               and not any(seg[t] == 3 and seg[t+1] == 3 for t in range(L-1)):
                out = F_seg(seg)[c - i]; break
        tab.append(out)
    return Rr, tab
for L in (1, 2, 3):
    Rr, tab = GL_table(L)
    print(f"G_{L} (radius {Rr}) commutes with F: {check_commutes_np_R(tab, F_loc, Rr)}; entries changed: {sum(1 for e,s in enumerate(tab) if s != index_word(e, 2*Rr+1)[Rr])}")
