# Steps 2-5 of the proof of kink rigidity for 94, end to end.  P = bare kink K^+_a (front 0).  W = a word w (tails 0 -> 1)
# after a front 0^b or 0^(b+1), a plateau of 1s and a bare descending kink at c (c != a mod 2).
# Lists the words w for which the outgoing Phi^+ part of (P, W) is P at both fronts; the theorem predicts the empty word only.
import sim, numpy as np, itertools
T=110; L=2*T+20; a=8; c=61
def ok94(s): return all(not(s[i]==0 and s[i+1]==1 and s[i+2]==0) for i in range(len(s)-2))
cores=[()]+[w for n in range(1,13) for w in itertools.product((0,1),repeat=n) if w[0]==1 and w[-1]==0 and ok94((0,0)+w+(1,1))]
print('cores',len(cores))
PI=np.array([0 if q<a else 1 for q in range(L)],np.uint8)
good=[]
for st in range(0,len(cores),2000):
    ch=cores[st:st+2000]; res=[]
    for b in (12,13):
        XI=[]
        for w in ch:
            s=[0]*b+list(w); s=s+[1]*(c-len(s))+[0]*(L-c); XI.append(np.array(s,np.uint8))
        cur,off,f=sim.run(94,[PI]*len(ch),XI,T,(0,0)); G,O=sim.outputs(cur,off,T,80,80); res.append(O)
    for k,w in enumerate(ch):
        if all((res[e][k][:80]==PI[:80]).all() for e in (0,1)): good.append(w)
print('D0 cores passing both shifts:',good[:10],len(good))
