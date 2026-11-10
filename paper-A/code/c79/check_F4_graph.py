# Confirmation (not used in the proof): (F^93_4) directly, by directed reachability in Gamma_4.
# Usage: python3 check_F4_graph.py k ...   (k=4: the five sets are strongly connected and closed; k=3: fails).
import numpy as np, sys, time
from graphs import build
from engine import vcell
def reach(n,src,dst,start):
    seen=np.zeros(n,bool); seen[start]=True
    while True:
        m=seen[src]; new=np.zeros(n,bool); new[dst[m]]=True; new&=~seen
        if not new.any(): return seen
        seen|=new
def run(k,umin=None,umax=None):
    umin=-k if umin is None else umin; umax=k if umax is None else umax
    valid,edges,Z,allc,K=build(k,umin,umax,ztr=[0])
    n=6*K
    src=np.concatenate([e[2] for e in edges]); dst=np.concatenate([e[3] for e in edges])
    # claimed class key: exterior cell; if exterior==0, tip==3 flag
    ext=np.tile(allc[:,0],6); tip=np.tile(allc[:,1],6)
    key=np.where(ext==0,4+(tip==3).astype(int),ext)
    keys=sorted(set(key[valid].tolist()))
    ok=True; info=[]
    for kk in keys:
        members=np.nonzero(valid&(key==kk))[0]
        s=members[0]
        fw=reach(n,src,dst,s); bw=reach(n,dst,src,s)
        cls=fw&bw
        # SCC of s must equal claimed set; and nothing outside reachable (fw subset claimed)
        eq=bool((cls==(valid&(key==kk))).all()) and bool((fw==(valid&(key==kk))).all())
        ok&=eq; info.append((kk,len(members),eq))
    return ok,info,int(valid.sum()),len(src)
if __name__=='__main__':
    for k in map(int,sys.argv[1:]):
        t0=time.time(); ok,info,nv,ne=run(k)
        print('F_%d'%k,'states',nv,'edges',ne,'claimed form holds (directed, all tails):',ok,info,flush=True)
        if k==4: assert ok and nv==4608 and ne==13824
        if k==3: assert sum(not e for _,_,e in info)==4
