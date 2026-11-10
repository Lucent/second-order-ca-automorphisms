# Exhaustive test of the one-exact-ray lemma (rule 94), front-0, front-1 and mirror versions, by direct simulation.
import sim, numpy as np, itertools, sys
NMAX=int(sys.argv[1]); T=70; L=2*T+20
def ok94(s): return all(not(s[i]==0 and s[i+1]==1 and s[i+2]==0) for i in range(len(s)-2))
def seqs(front,nfront,nmax,first=None):
    out=[]
    for n in range(0,nmax+1):
        for w in itertools.product((0,1),repeat=n):
            if first is not None and (n==0 or w[0]!=first): continue
            for tail in (0,1):
                s=[front]*nfront+list(w)
                full=s+[tail]*(L-len(s))
                if ok94(full): out.append(full)
    # dedupe
    u={tuple(x) for x in out}; return [list(x) for x in u]
def k10(s,par):  # 1->0 kinks x with x == par mod 2
    return [x for x in range(1,len(s)) if s[x-1]==1 and s[x]==0 and (x-par)%2==0]
def test(front,mirror):
    bad=0; tot=0; badloc=0
    for a1 in (5,6):
        As=seqs(front,a1,NMAX,first=1-front)   # A=front before a1, A(a1)=1-front
        Bs=[]
        for b0 in (2,3):
            Bs+=seqs(front,b0,NMAX)
        Bs=[list(x) for x in {tuple(x) for x in Bs}]
        pairs=[(a,b) for a in As for b in Bs]
        for s0 in range(0,len(pairs),4000):
            ch=pairs[s0:s0+4000]
            if not mirror: PI=[np.array(a,np.uint8) for a,b in ch]; XI=[np.array(b,np.uint8) for a,b in ch]
            else: PI=[np.array(b,np.uint8) for a,b in ch]; XI=[np.array(a,np.uint8) for a,b in ch]
            cur,off,f=sim.run(94,PI,XI,T,(front,front))
            for c,(a,b) in enumerate(ch):
                if not mirror: ray=[x for x in f[c] if x[0]==a1-1]; pos=lambda x:x[1]
                else: ray=[x for x in f[c] if x[1]==a1-1]; pos=lambda x:x[0]
                ks=k10(b,a1)
                pred=1 if ks else 0
                tot+=1
                if len(ray)!=pred: bad+=1
                elif pred and pos(ray[0])!=ks[0]-1: badloc+=1
    print('front',front,'mirror',mirror,'cases',tot,'count mismatches',bad,'location mismatches',badloc,flush=True)
for front in (0,1):
    for mirror in (False,True):
        test(front,mirror)
