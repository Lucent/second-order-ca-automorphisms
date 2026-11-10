from strip import tree
import sys
K=int(sys.argv[1])
levels=tree(0,K,verbose=False)
cyc=levels[-1]
print('width',K,'00-front cycles',len(cyc))
phim=0; sec=0; last=0; other=0
def front11(s): return s.endswith('011') and '11' not in s[:-1] and '000' not in s
for c in cyc:
    L=len(c)
    const=[all(c[t][j]==c[0][j] for t in range(L)) for j in range(K)]
    if all(const):
        s='10'+''.join(str(v) for v in c[0])
        if '11' not in s and '000' not in s: phim+=1
        elif front11(s): last+=1
        else: other+=1
        continue
    j1=const.index(False)
    s='10'+''.join(str(c[0][j]) for j in range(j1))
    if front11(s): sec+=1
    else: other+=1
print('Phi- words',phim,'shifted 11-fronts',sec,'shifted 11-fronts with the final 11 on the last two lines',last,'other',other)
