import sys,time,json,argparse
from moves_engine import *
ap=argparse.ArgumentParser()
ap.add_argument('rule',type=int); ap.add_argument('--bx',default='0'); ap.add_argument('--by',default='0')
ap.add_argument('--speeds',default='1/1'); ap.add_argument('--kmin',type=int,default=1); ap.add_argument('--kmax',type=int,default=5)
ap.add_argument('--tlo',type=int,default=4); ap.add_argument('--thi',type=int,default=2); ap.add_argument('--S',type=int,default=0)
ap.add_argument('--L0',type=int,default=6); ap.add_argument('--maxk',type=int,default=6); ap.add_argument('--out',default=''); ap.add_argument('--minclean',type=float,default=0.0); ap.add_argument('--xext',type=int,default=0); ap.add_argument('--maxscc',type=int,default=10**9)
A=ap.parse_args()
speeds=[tuple(map(int,x.split('/'))) for x in A.speeds.split(',')]
bg=(tuple(int(c) for c in A.bx),tuple(int(c) for c in A.by))
E=Eng(A.rule,bg,speeds,maxk_per_state=A.maxk,maxscc=A.maxscc)
import numpy as np
import moves_engine as ufeng
ufeng.XEXT=bool(A.xext)
S=A.S or 10+8*max(b for a,b in speeds)
# prefilter kinks: keep those with some clean bounce at k=3,t0 in {0,-1,-2}
def prefilter(E,S):
    k=3; K4=4**k; N=4*K4; P=2*S+40; r=P+k; Wd=r+1+2*S+60
    rng=np.random.default_rng(5); n=np.arange(N); sidx=n%K4; ext=(n//K4)%4
    for st0 in range(len(E.ALL)):
        X0=np.zeros((Wd,N//64),np.uint64);Y0=np.zeros((Wd,N//64),np.uint64)
        for c in range(P): X0[c]=pack(rng.integers(0,2,N));Y0[c]=pack(rng.integers(0,2,N))
        X0[P]=pack(ext>>1);Y0[P]=pack(ext&1)
        for j in range(k): X0[P+1+j]=pack((sidx>>(2*j+1))&1);Y0[P+1+j]=pack((sidx>>(2*j))&1)
        E.fill_state(X0,Y0,r+1,Wd,E.ALL[st0])
        keep=[]
        for kk in E.lib[st0]:
            tot=0
            for t0 in (0,-1,-2):
                X,Y,valid,ns=E.move(X0,Y0,st0,r,0,kk,S,A.L0,gap=kk['gap']); tot+=unpack(valid,N).sum()
            kk['clean_frac']=round(tot/(3*N),3)
            if tot>0 and kk['clean_frac']>=A.minclean: keep.append(kk)
        E.lib[st0]=keep
prefilter(E,S)
print('rule',A.rule,'bg',bg,'states',E.ALL,'rho_ok',E.rho_ok,'S',S)
for i,L in E.lib.items(): print('  lib state',i,[(x['a'],x['b'],x['path'][:10],x['cyc'],'gap',x['gap'],x['clean_frac']) for x in L])
res=[]
for k in range(A.kmin,A.kmax+1):
    t=time.time()
    lab,st,r=run_level(E,k,-k-A.tlo-(S if False else 0),k+A.thi,S=S,L0=A.L0)
    an=analyse(E,k,lab,r)
    used=sorted(set((x[0],x[1]) for x in st.get('used_moves',[])))
    print('k',k,'edges',st['edges'],'ctxdis',st['ctxdis'],'used t0 (bar0)',[x[1] for x in used if x[0]==0],'(bar1)',[x[1] for x in used if x[0]==1],'%.1fs'%(time.time()-t))
    for ex,v in an.items(): print('    ext',ex,v)
    sys.stdout.flush()
    res.append(dict(k=k,stats={kk:vv for kk,vv in st.items() if kk not in('used_moves','dropped_moves')},used=used,an=an))
if A.out: json.dump(res,open(A.out,'w'))
