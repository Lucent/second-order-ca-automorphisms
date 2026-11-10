# d'Alembert bookkeeping: for E0 (front beta) vs arbitrary B (front beta, any back tail),
# o(q) xor pi(q) xor [B(inf) xor kappa] == parity of firings on rays q'<q, q' != q mod 2 ; check all q.
import sim, numpy as np, random
random.seed(4); T=90; L=2*T+20
def ok94(s): return all(not(s[i]==0 and s[i+1]==1 and s[i+2]==0) for i in range(len(s)-2))
def rnd(beta):
    while True:
        n=random.randint(0,14); w=[random.randint(0,1) for _ in range(n)]; tail=random.randint(0,1)
        s=[beta]*random.randint(2,5)+w; s=s+[tail]*(L-len(s))
        if ok94(s): return s
bad=0; N=0
for beta in (0,1):
    cs=[(rnd(beta),rnd(beta)) for _ in range(400)]
    cur,off,f=sim.run(94,[np.array(a,np.uint8) for a,b in cs],[np.array(b,np.uint8) for a,b in cs],T,(beta,beta))
    G,O=sim.outputs(cur,off,T,50,50)
    for c,(a,b) in enumerate(cs):
        N+=1
        for q in range(2,45):
            n=sum(1 for (qq,pp) in f[c] if qq<q and (qq-q)%2)%2
            if (O[c][q]^a[q]^b[-1]^beta)!=n: bad+=1; break
print('cases',N,'bad',bad)
