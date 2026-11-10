"""Rules 127R, 179R, 223R, Proposition 6.19: the family W generated from one word under F and F^-1 in all
contexts is a single orbit of context-free words; the involution iota (half a period) changes only
position 2; (H1) in all contexts; on all words of length up to 10 that begin with a word of W, every flip
preserves the occurrences and every two flips commute (H2); F^P fixes each word in random contexts.
"""
import numpy as np, itertools
from hw_ca import *
ex={127:'03112',179:'01110',223:'00122'}
k=5
def W_of(rule,w0):
    Ft,Fi=tables(rule); W={w0}; fr=[w0]
    while fr:
        w=fr.pop()
        for T in (Ft,Fi):
            for l,r in itertools.product(range(4),repeat=2):
                u=step_word(T,(l,)+w+(r,))
                if u not in W: W.add(u); fr.append(u)
    return W,Ft,Fi
res={}
for rule,s in ex.items():
    w0=tuple(int(ch) for ch in s); W,Ft,Fi=W_of(rule,w0)
    # context freeness: image unique
    cf=all(len({step_word(T,(l,)+w+(r,)) for l in range(4) for r in range(4)})==1 for w in W for T in (Ft,Fi))
    nxt={w:step_word(Ft,(0,)+w+(0,)) for w in W}
    orb=[w0]
    while nxt[orb[-1]]!=w0: orb.append(nxt[orb[-1]])
    P=len(orb); h=P//2
    iot={orb[i]:orb[(i+h)%P] for i in range(P)}
    diffpos={tuple(j for j in range(k) if w[j]!=iot[w][j]) for w in W}
    # time-periodic: F^P fixes each word in every context: follows from cf + single orbit; also check directly with random contexts
    rng=np.random.default_rng(rule); tpok=True
    for w in W:
        for t in range(300):
            c=rng.integers(0,4,2*P+k+4); st=P+2; c[st:st+k]=w
            x=c.copy()
            # use a big ring with random context to emulate infinite context for P steps
            big=rng.integers(0,4,400); big[100:100+len(c)]=c
            for _ in range(P): big=step(Ft,big)
            if tuple(big[100+st:100+st+k])!=w: tpok=False
    # H1 exhaustive: flip commutes in every context (7-cell windows enough since iota changes only interior positions)
    h1=True
    for w in W:
        for l,r in itertools.product(range(4),repeat=2):
            for T in (Ft,Fi):
                a=step_word(T,(l,)+w+(r,)); b=step_word(T,(l,)+iot[w]+(r,))
                if b!=iot[a]: h1=False
    # need also cells just outside the window: they read position 0 / k-1, check iota fixes those
    edge=all(w[0]==iot[w][0] and w[-1]==iot[w][-1] for w in W)
    # H2 exhaustive over words of length up to 2k+? (pairs at offset 1..k-1, plus non-overlapping)
    def occ(u): return [i for i in range(len(u)-k+1) if tuple(u[i:i+k]) in W]
    def fl(u,i): u=list(u); u[i:i+k]=iot[tuple(u[i:i+k])]; return tuple(u)
    h2=0; offs=set()
    for n in range(k+1,2*k+1):
        # enumerate words with an occurrence at 0 and another occurrence: build from W prefix
        for w in W:
            for rest in itertools.product(range(4),repeat=n-k):
                u=w+rest; o=occ(u)
                if len(o)<2: continue
                for i in o:
                    if occ(fl(u,i))!=o: h2+=1
                    for j in o:
                        if j>i:
                            offs.add(j-i)
                            if fl(fl(u,i),j)!=fl(fl(u,j),i): h2+=1
    res[rule]=(W,orb,iot,Ft,Fi)
    print(rule,'|W|=',len(W),'single orbit P=',P,P==len(W),'context-free F,F^-1:',cf,'F^P fixes (random contexts):',tpok,
          'iota diff positions',diffpos,'H1:',h1,'edges agree:',edge,'H2 failures',h2,'overlap offsets seen',sorted(offs))
    print('   orbit',[''.join(map(str,w)) for w in orb])
