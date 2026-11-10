import numpy as np, random, sys
from core import *
def ext(cyc,beta,k):
    """cyc: list of tuples (Z_0..Z_{k-1}) at t=0..L-1 (periodic). returns list of (kind, newcyc list)"""
    L=len(cyc)
    def line(t,j):
        if j==-1: return beta
        if j==-2: return 1-beta
        return cyc[t%L][j]
    def run(z,reps):
        out=[]
        for t in range(reps*L):
            out.append(z)
            z=z^line(t-1,k-2)^line(t,k-2)^N(line(t,k-2),line(t,k-1),z)
        return out,z
    ends={z:run(z,1)[1] for z in (0,1)}
    res=[]
    if ends[0]==ends[1]:
        z=ends[0]; seq,_=run(z,1)   # start on the cycle: after one period, value at t=0 is ends
        # check consistency: run(z)[1]==z
        assert run(z,1)[1]==z
        res.append(('slaved',[cyc[t]+(seq[t],) for t in range(L)]))
    elif ends[0]==0:
        for z in (0,1):
            seq,_=run(z,1); res.append(('free',[cyc[t]+(seq[t],) for t in range(L)]))
    else:
        seq,e=run(0,2); assert e==0
        res.append(('double',[cyc[t%L]+(seq[t],) for t in range(2*L)]))
    return res
def tree(beta,K,verbose=True):
    cyc=[[(beta,)]]
    out=[]
    for k in range(1,K):
        new=[];kinds={}
        for c in cyc:
            r=ext(c,beta,k)
            kd=r[0][0]; kinds[kd]=kinds.get(kd,0)+1
            new+= [x[1] for x in r]
        cyc=new
        if verbose: print('line',k,'branches',len(cyc),kinds,'periods',sorted(set(len(c) for c in cyc)))
        out.append(cyc)
    return out
if __name__=='__main__':
    tree(int(sys.argv[1]),int(sys.argv[2]))
