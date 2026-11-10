# Brute-force scattering closure with the local property (C0), from the bare points only (rule 94R).
# Universe: sequences lam^inf core tail^inf with |core|<=CAP, stored as (lam, core); a stored sequence is
# used both as a point of Phi^- and (by rev) of Phi^+.  Checks that every object with |core|<=N is reached.
# usage: python3 closure_bf.py L N CAP PCAP   (window length L; probes are the sequences with |core|<=PCAP)
from lib import scatter, core as corefn, inphi
import sys, itertools
L=int(sys.argv[1]); N=int(sys.argv[2]); CAP=int(sys.argv[3]); PCAP=int(sys.argv[4])
def objs(n):
    for lam in '01':
        for bits in itertools.product('01',repeat=n):
            c=''.join(bits)
            if c[0]==lam: continue
            if n>1 and c[-1]==c[-2]: continue
            if not inphi(lam*3+c+c[-1]*3): continue
            yield (lam,c)
allo=[('0',''),('1','')]+[o for n in range(1,CAP+1) for o in objs(n)]
def words(o):
    lam,c=o; s=lam*L+c+(c[-1] if c else lam)*L
    return {s[i:i+L] for i in range(len(s)-L+1)}
def mu(o):
    lam,c=o
    if not c: return o
    return (c[-1],c[::-1][1:]+lam)
def S(P,X,par):
    # P=(lam,core) as probe, X=(lam,core) as object, same lam; relative placement parity par
    lam=P[0]
    pw=lam*(2+par)+P[1] if P[1] else lam*2
    xw=lam*4+X[1] if X[1] else lam*2
    g,o,F=scatter(pw,xw)
    return corefn(g),corefn(o)
known={('0',''),('1',''),('0','1'),('1','0')}
W=set()
for o in known: W|=words(o)
rnd=0
cache={}
while True:
    rnd+=1; new=set()
    for o in allo:                       # (C0)
        if o not in known and words(o)<=W: new.add(o)
    for o in list(known):                # (C1): mu (rev is built in, shifts are built in)
        m=mu(o)
        if m not in known: new.add(m)
    probes=[p for p in allo if len(p[1])<=PCAP]
    for P in probes:
        for X in allo:
            if P[0]!=X[0]: continue
            if not P[1] or not X[1]: continue
            for par in (0,1):
                key=(P,X,par)
                if key not in cache: cache[key]=S(P,X,par)
                G,O=cache[key]
                if P in known and X in known:                 # (C2)
                    for z in (G,O):
                        if z not in known and len(z[1])<=CAP: new.add(z)
                elif G in known and O in known:               # (C3)
                    for z in (P,X):
                        if z not in known: new.add(z)
    if not new: break
    known|=new
    for o in new: W|=words(o)
    small=sum(1 for o in allo if len(o[1])<=N)
    ksmall=sum(1 for o in allo if len(o[1])<=N and o in known)
    print('round',rnd,'known',len(known),'of',len(allo),'; core<=',N,':',ksmall,'/',small,flush=True)
    if ksmall==small: break
missing=[o for o in allo if len(o[1])<=N and o not in known]
print('missing',missing[:20])
