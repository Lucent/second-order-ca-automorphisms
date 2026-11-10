# Rule 95: (1) f_95 does not read its middle argument; (2) the two parts of the row lemma, by direct simulation;
# (3) all pairs of cores that satisfy both relations (T) and (Delta) of the proof of kink rigidity for 95.
# usage: t95b.py K
import sim, numpy as np, itertools, sys
t=sim.ftab(95); print('f95 = not(a and c):', all(t[4*a+2*b+c]==1-(a&c) for a in (0,1) for b in (0,1) for c in (0,1)))
K=int(sys.argv[1]); T=3*K+40; L=2*T+20
def ok95(s): return all(not(s[i]==0 and s[i+2]==0) for i in range(len(s)-2))
cores=[c for n in range(1,K+1) for c in itertools.product((0,1),repeat=n) if c[0]==0 and c[-1]==0 and ok95((1,1)+c+(1,1))]
print('cores',len(cores))
def place(c,s): a=np.ones(L,np.uint8); a[s:s+len(c)]=c; return a
# first/second row lemma on Y level: E with first zero (on sublattice s) at a1; D arbitrary.
# prediction: o(a1)=0 (defect survives) iff D has no zero x with x==a1 mod 2;
#             o(a1+2)=0 iff [E(a1+2)=1 and D has exactly one zero ==a1]  -- second-row (if E(a1+2)=0, row a1+2 starts its own carry; check only E(a1+2)=1)
bad1=bad2=n2=0; tot=0
pairs=[(ce,cd,s) for ce in cores for cd in cores for s in (10,11)]
sols=[]
for st in range(0,len(pairs),3000):
    ch=pairs[st:st+3000]
    PI=[place(ce,10) for ce,cd,s in ch]; XI=[place(cd,s) for ce,cd,s in ch]
    cur,off,f=sim.run(95,PI,XI,T,(1,1)); G,O=sim.outputs(cur,off,T,L//2,L//2)
    for c,(ce,cd,s) in enumerate(ch):
        A=PI[c]; B=XI[c]
        for par in (0,1):
            za=[x for x in range(L) if A[x]==0 and x%2==par]
            if not za: continue
            a1=za[0]; zd=[x for x in range(L) if B[x]==0 and x%2==par]
            tot+=1
            if (O[c][a1]==0)!=(len(zd)==0): bad1+=1
            if A[a1+2]==1:
                n2+=1
                if (O[c][a1+2]==0)!=(len(zd)==1): bad2+=1
        # relations: s=11 odd placement (B start parity differs from A start 10)
    # collect relations
    for c,(ce,cd,s) in enumerate(ch):
        pass
print('first-row tests',tot,'mismatch',bad1,'second-row tests',n2,'mismatch',bad2)
# full relation search: transparent at one placement parity, both delayed by 2 at the other
def S(ce,cd,s):
    cur,off,f=sim.run(95,[place(ce,10)],[place(cd,s)],T,(1,1)); G,O=sim.outputs(cur,off,T,L//2,L//2); return G[0],O[0]
found=[]
for st in range(0,len(cores)**2,3000):
    ch=[(ce,cd) for ce in cores for cd in cores][st:st+3000]
    res={}
    for s in (10,11,12,13):
        PI=[place(ce,10) for ce,cd in ch]; XI=[place(cd,s) for ce,cd in ch]
        cur,off,f=sim.run(95,PI,XI,T,(1,1)); res[s]=sim.outputs(cur,off,T,L//2,L//2)
    M=L//2-T//2
    for c,(ce,cd) in enumerate(ch):
        for so,se in ((11,10),(10,11)):  # (odd-placement, even-placement) shifts of D relative to E
            Go,Oo=res[so][0][c],res[so][1][c]; Ge,Oe=res[se][0][c],res[se][1][c]
            trans = (Go[:M]==place(cd,so)[:M]).all() and (Oo[:M]==place(ce,10)[:M]).all()
            dly = (Ge[:M]==place(cd,se+2)[:M]).all() and (Oe[:M]==place(ce,12)[:M]).all()
            if trans and dly: found.append((ce,cd,so,se))
print('solutions of both relations:',found)
