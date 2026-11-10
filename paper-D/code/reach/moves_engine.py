# Exact cone-base orbits of edge moves built from REAL probes (bit-packed over the batch).
import numpy as np, ctypes, os, sys, json, time, random
from common import *
from sim import Fst, bg_orbit
from sft import kinks
from bgstates import shift, rho
_lib=ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),'uf.so'))
def ufind(n,u,v):
    p=np.arange(n,dtype=np.int64); u=np.ascontiguousarray(u,np.int64); v=np.ascontiguousarray(v,np.int64)
    P=ctypes.POINTER(ctypes.c_int64)
    _lib.uf(ctypes.c_int64(n),ctypes.c_int64(len(u)),u.ctypes.data_as(P),v.ctypes.data_as(P),p.ctypes.data_as(P))
    return p
ONES=np.uint64(0xFFFFFFFFFFFFFFFF)
class Eng:
    def __init__(s,rule,bg0,speeds,maxk_per_state=6,kfilter=None,maxscc=10**9):
        s.rule=rule; s.mins=[t for t in range(8) if (rule>>t)&1]
        orb=bg_orbit(rule,bg0); s.T=len(orb); s.p=len(bg0[0])
        # states as (ph,d)
        s.ALL=[]; s.key={}
        for ph,o in enumerate(orb):
            for d in range(s.p):
                st=shift(o,d)
                if st in s.key: continue
                s.key[st]=len(s.ALL); s.ALL.append(st)
        s.orb=orb
        s.Fidx=[s.key[Fst(rule,st)] for st in s.ALL]
        s.Bidx=[s.Fidx.index(i) for i in range(len(s.ALL))]
        s.rhoidx=[s.key.get(rho(st),-1) for st in s.ALL]
        s.rho_ok=all(i>=0 for i in s.rhoidx)
        # probe library per state
        s.lib={}
        for i,st in enumerate(s.ALL):
            L=[]
            for (a,b) in speeds:
                kk=kinks(rule,a,b,st,maxpaths=MAXPATHS)
                if kk is None: continue
                for k in kk[1]:
                    path=[int(c) for c in k['path']]
                    cyc=None if k['cycle'] is None else [int(c) for c in k['cycle']]
                    if kfilter and not kfilter(a,b,path,cyc): continue
                    if maxscc<10**9 and (k['scc_size']>maxscc or (cyc is not None and not k['scc_simple'])): continue
                    if any(x['path']==path and x['cyc']==cyc for x in L): continue
                    kd=dict(a=a,b=b,path=path,cyc=cyc,j=k['Bphase'],gend=k.get('glider_end'))
                    from legit import legit_gap
                    g=legit_gap(rule,kd,st,s.p)
                    if g is None: continue
                    kd['gap']=g; L.append(kd)
            s.lib[i]=L[:maxk_per_state]
    def advance(s,i,n):
        for _ in range(abs(n)): i=s.Fidx[i] if n>0 else s.Bidx[i]
        return i
    def f(s,a,b,c):
        o=np.zeros_like(b)
        for t in s.mins:
            o|=(a if t&4 else ~a)&(b if t&2 else ~b)&(c if t&1 else ~c)
        return o
    def go(s,X,Y,n):
        for _ in range(abs(n)):
            if n>0:
                X,Y=Y,X^s.f(np.roll(Y,1,0),Y,np.roll(Y,-1,0))
            else:
                X,Y=Y^s.f(np.roll(X,1,0),X,np.roll(X,-1,0)),X
        return X,Y
    def fill_state(s,X,Y,lo,hi,st):
        xs,ys=st
        for i in range(lo,hi):
            X[i]=ONES if xs[i%s.p] else 0; Y[i]=ONES if ys[i%s.p] else 0
    def place(s,X,Y,kk,st_idx,start_min):
        p=s.p; L=2*kk['b']; j=kk['j']
        st0=start_min
        while (st0-(j+L))%p: st0+=1
        i=st0; Wd=X.shape[0]
        for c in kk['path']:
            if i>=Wd: break
            X[i]=ONES if c>>1 else 0; Y[i]=ONES if c&1 else 0; i+=1
        if kk['cyc'] is None:
            xs,ys=s.ALL[st_idx]; cells=[2*xs[q]+ys[q] for q in range(p)]; o=(kk['gend']+L)%p
            while i<Wd:
                c=cells[o%p];X[i]=ONES if c>>1 else 0; Y[i]=ONES if c&1 else 0; i+=1;o+=1
        else:
            jj=0
            while i<Wd:
                c=kk['cyc'][jj%len(kk['cyc'])];X[i]=ONES if c>>1 else 0; Y[i]=ONES if c&1 else 0; i+=1;jj+=1
    def move(s,X0,Y0,st0,r,t0,kk,S,L0,gap=2):
        """apply probe kk glued at time t0 at r+gap. X0,Y0 time-0 state (tail st0). returns X,Y,valid,newst(array of state idx)"""
        X,Y=s.go(X0.copy(),Y0.copy(),t0)
        stg=s.advance(st0,t0)
        X=X.copy();Y=Y.copy(); s.fill_state(X,Y,r+1,X.shape[0],s.ALL[stg])
        s.place(X,Y,kk,stg,r+gap)
        X,Y=s.go(X,Y,S)
        valid=np.zeros(X.shape[1],np.uint64); newst=[]
        Xc=X.copy();Yc=Y.copy()
        for i,st in enumerate(s.ALL):
            xs,ys=st; m=ONES*np.ones(X.shape[1],np.uint64)
            for c in range(r+1,r+1+L0):
                m&=(X[c] if xs[c%s.p] else ~X[c])&(Y[c] if ys[c%s.p] else ~Y[c])
            edgeB=(X[r] if xs[r%s.p] else ~X[r])&(Y[r] if ys[r%s.p] else ~Y[r])
            m&=~edgeB; m&=~valid
            if not m.any(): continue
            newst.append((i,m)); valid|=m
            for c in range(r+1,X.shape[0]):
                bx=ONES if xs[c%s.p] else np.uint64(0); by=ONES if ys[c%s.p] else np.uint64(0)
                Xc[c]=(Xc[c]&~m)|(bx&m); Yc[c]=(Yc[c]&~m)|(by&m)
        X,Y=s.go(Xc,Yc,-(t0+S))
        ns=[(s.advance(i,-(t0+S)),m) for i,m in newst]
        return X,Y,valid,ns
def unpack(a,N): return np.unpackbits(a.view(np.uint8),bitorder='little')[:N]
def pack(bits): return np.packbits(bits.astype(np.uint8),bitorder='little').view(np.uint64)
STRICT=True
MAXPATHS=400
XEXT=False
RECORD=None
def run_level(E,k,tmin,tmax,S,L0=6,seed=1,bars=(0,1),log=print,NCTX=4):
    rule=E.rule; nS=len(E.ALL); K4=4**k
    NSUB=4*K4; N=max(64,NCTX*NSUB); N=((N+63)//64)*64
    Tmax=max(abs(tmin),abs(tmax))
    ntot=Tmax+S+Tmax+S
    P=ntot+4; r=P+k; Rr=ntot+L0+40; Wd=r+1+Rr
    rng=np.random.default_rng(seed)
    n=np.arange(N); sidx=n%K4; ext=(n//K4)%4
    nodes_u=[];nodes_v=[]
    stats=dict(moves=0,edges=0,ctxdis=0,leak=0,unclean=0)
    def nid(st,e,sw): return (st*4+e)*K4+sw
    for st0 in range(nS):
        X0=np.zeros((Wd,N//64),np.uint64);Y0=np.zeros((Wd,N//64),np.uint64)
        for c in range(P): X0[c]=pack(rng.integers(0,2,N));Y0[c]=pack(rng.integers(0,2,N))
        X0[P]=pack(ext>>1);Y0[P]=pack(ext&1)
        for j in range(k):
            X0[P+1+j]=pack((sidx>>(2*j+1))&1);Y0[P+1+j]=pack((sidx>>(2*j))&1)
        E.fill_state(X0,Y0,r+1,Wd,E.ALL[st0])
        xs,ys=E.ALL[st0]; ebits=(2*xs[r%E.p]+ys[r%E.p])
        edge_cell=(sidx>>(2*(k-1)))&3
        okstart=(edge_cell!=ebits)
        for bar in bars:
            if bar and not E.rho_ok: continue
            if bar: Xs,Ys=Y0,X0; sts=E.rhoidx[st0]
            else: Xs,Ys=X0,Y0; sts=st0
            for t0 in range(tmin,tmax+1):
                stg=E.advance(sts,t0)
                for kk in E.lib[stg]:
                    X,Y,valid,ns=E.move(Xs,Ys,sts,r,t0,kk,S,L0,gap=kk['gap'])
                    if bar: X,Y=Y,X; ns=[(E.rhoidx[i],m) for i,m in ns]
                    stats['moves']+=1
                    # unchanged outside window (exact region ntot+1..P)
                    ch=np.zeros(N//64,np.uint64)
                    for c in range(ntot+1,P): ch|=(X[c]^X0[c])|(Y[c]^Y0[c])
                    ch|=(Y[P]^Y0[P])
                    if not XEXT: ch|=(X[P]^X0[P])
                    nex=(unpack(X[P],N).astype(np.int64)<<1)|unpack(Y[P],N).astype(np.int64)
                    v=unpack(valid,N).astype(bool); chb=unpack(ch,N).astype(bool)
                    newst=np.full(N,-1)
                    for i,m in ns: newst[unpack(m,N).astype(bool)]=i
                    sw=np.zeros(N,np.int64)
                    for j in range(k):
                        sw|=unpack(X[P+1+j],N).astype(np.int64)<<(2*j+1); sw|=unpack(Y[P+1+j],N).astype(np.int64)<<(2*j)
                    ok=v&~chb&okstart
                    allok=np.ones(NSUB,bool); agree=np.ones(NSUB,bool); anyok=np.zeros(NSUB,bool)
                    for cc in range(NCTX):
                        sl=slice(cc*NSUB,(cc+1)*NSUB); allok&=ok[sl]; anyok|=ok[sl]
                        agree&=(sw[sl]==sw[:NSUB])&(newst[sl]==newst[:NSUB])&(nex[sl]==nex[:NSUB])
                    good=allok&agree
                    if STRICT:
                        lk=np.zeros(NSUB,bool)
                        for cc in range(NCTX):
                            sl=slice(cc*NSUB,(cc+1)*NSUB); lk|=(v[sl]&chb[sl])
                        lk_s=lk.reshape(4,K4).any(axis=0)   # leak for some ext at this window
                        good&=~np.tile(lk_s,4)
                    nd=int((allok&~agree).sum())+int((anyok&~allok).sum()); stats['ctxdis']+=nd
                    if nd>0:
                        stats.setdefault('dropped_moves',[]).append((bar,t0)); continue
                    stats.setdefault('used_moves',[]).append((bar,t0,int(good.sum()))); stats['leak']+=int((okstart[:NSUB]&v[:NSUB]&chb[:NSUB]).sum())
                    stats['unclean']+=int((okstart[:NSUB]&~v[:NSUB]).sum())
                    if RECORD is not None:
                        # clean-in-all-contexts & context-free result (leaks allowed: record leak flag)
                        cl=np.ones(NSUB,bool); ag=np.ones(NSUB,bool); lk=np.zeros(NSUB,bool)
                        for cc in range(NCTX):
                            sl=slice(cc*NSUB,(cc+1)*NSUB); cl&=(v[sl]&okstart[sl]); lk|=chb[sl]
                            ag&=(sw[sl]==sw[:NSUB])&(newst[sl]==newst[:NSUB])
                        tgt=np.where(cl&ag&~lk,nid(newst[:NSUB],nex[:NSUB],sw[:NSUB]),np.where(cl&~lk,-3,np.where(cl,-2,-1)))
                        RECORD.setdefault((bar,t0,str(kk['path'])),{})[st0]=tgt
                    idx=np.nonzero(good)[0]
                    u=nid(st0,ext[idx],sidx[idx]); w=nid(newst[idx],nex[idx],sw[idx])
                    sel=u!=w
                    nodes_u.append(u[sel]);nodes_v.append(w[sel]); stats['edges']+=int(sel.sum())
    Nn=nS*4*K4
    lab=ufind(Nn,np.concatenate(nodes_u) if nodes_u else np.zeros(0,np.int64),np.concatenate(nodes_v) if nodes_v else np.zeros(0,np.int64))
    return lab,stats,r
def analyse(E,k,lab,r):
    nS=len(E.ALL); K4=4**k; f=E.rule; p=E.p
    Pa={(b,c) for b in (0,1) for c in (0,1) if fv(f,0,b,c)!=fv(f,1,b,c)}
    nodes=np.arange(nS*4*K4); st=nodes//(4*K4); e=(nodes//K4)%4; s=nodes%K4
    Bcell=np.array([2*E.ALL[i][0][r%p]+E.ALL[i][1][r%p] for i in range(nS)])
    valid=((s>>(2*(k-1)))&3)!=Bcell[st]
    yt=s&1
    if k>=2: yt1=(s>>2)&1
    else: yt1=np.array([E.ALL[i][1][(r+1)%p] for i in range(nS)])[st]
    read=np.zeros(len(nodes),bool)
    for (b,c) in Pa: read|=(yt==b)&(yt1==c)
    read&=valid
    res={}
    for ex in range(4):
        sel=valid&(e==ex)
        L=lab[sel]; u=np.unique(L)
        hit=np.unique(lab[read&(e==ex)])
        bad=np.setdiff1d(u,hit)
        nbadstates=int(np.isin(L,bad).sum())
        # orbit sizes
        _,cnt=np.unique(L,return_counts=True)
        res[ex]=dict(states=int(sel.sum()),orbits=len(u),noread_orbits=len(bad),noread_states=nbadstates,sizes=sorted(cnt.tolist())[:6] if len(cnt)<=64 else 'many')
    return res
