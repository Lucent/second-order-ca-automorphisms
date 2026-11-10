# Closure that forces only the value at the point S: window w (centre (0,0)) is forced if for some extension (D=2 cells/side, rows -1,0)
# the 5-point stencil centred at (1,0) [points N=(2,0), S=(0,0), W=(1,-1), C=(1,0), E=(1,1)] admits no G with
# G(S)=1 and G=0 at every point whose window is already forced.
import numpy as np, sys, itertools
from core import FT
def run(R,Z0,maxit=200,verbose=False):
    W=2*R+1; NW=4**W; L0=W+4
    w=np.arange(NW,dtype=np.int64)
    Xc=np.zeros((NW,W),np.uint8); Yc=np.zeros((NW,W),np.uint8)
    for k in range(W):
        Xc[:,k]=(w>>(W+W-1-k))&1; Yc[:,k]=(w>>(W-1-k))&1
    pw=(1<<np.arange(W-1,-1,-1)).astype(np.int64)
    def widx(xp,yp): return (xp.astype(np.int64)@pw<<W)|(yp.astype(np.int64)@pw)
    pre=[]
    for e in range(256):
        b=[(e>>j)&1 for j in range(8)]
        X=np.concatenate([np.tile(np.array(b[0:2],np.uint8),(NW,1)),Xc,np.tile(np.array(b[2:4],np.uint8),(NW,1))],1)
        Y=np.concatenate([np.tile(np.array(b[4:6],np.uint8),(NW,1)),Yc,np.tile(np.array(b[6:8],np.uint8),(NW,1))],1)
        # columns 0..L0-1 <-> positions -R-2..R+2
        Y1=X[:,1:-1]^FT[4*Y[:,:-2]+2*Y[:,1:-1]+Y[:,2:]]          # positions -R-1..R+1
        Y2=Y[:,2:-2]^FT[4*Y1[:,:-2]+2*Y1[:,1:-1]+Y1[:,2:]]       # positions -R..R
        # windows at (1,i): x=Y0, y=Y1, positions i-R..i+R
        def win1(i):
            a=i+1  # Y1 column of position i-R is (i-R)+(R+1)
            return widx(Y[:,a+1:a+1+W], Y1[:,a:a+W])
        iW=win1(-1); iC=win1(0); iE=win1(1); iN=widx(Y1[:,1:1+W],Y2)
        yW=Y1[:,R]; yC=Y1[:,R+1]; yE=Y1[:,R+2]
        pre.append((iN.astype(np.int32),iW.astype(np.int32),iC.astype(np.int32),iE.astype(np.int32),yW,yC,yE))
    Z=Z0.copy()
    for it in range(maxit):
        newZ=Z.copy()
        for (iN,iW,iC,iE,yW,yC,yE) in pre:
            fN,fW,fC,fE=Z[iN],Z[iW],Z[iC],Z[iE]
            base=FT[4*yW+2*yC+yE]
            sat=np.zeros(NW,bool)
            for gN,gW,gC,gE in itertools.product((0,1),repeat=4):
                ok=np.ones(NW,bool)
                if gN: ok&=~fN
                if gW: ok&=~fW
                if gC: ok&=~fC
                if gE: ok&=~fE
                rel=((gN^1)==(FT[4*(yW^gW)+2*(yC^gC)+(yE^gE)]^base))
                sat|=ok&rel
            newZ|=~sat
        if verbose: print(' iter',it,int(newZ.sum()),flush=True)
        if (newZ==Z).all(): break
        Z=newZ
    return NW-int(Z.sum()), it
if __name__=='__main__':
    from seeds import base, fronts
    R=int(sys.argv[1]); mode=sys.argv[2]
    W=2*R+1; B=base(R)
    if mode=='none': S=B
    elif mode=='file': 
        from seeds import load; S=load(sys.argv[3],R)
    else:
        betas=tuple(int(c) for c in sys.argv[3]); rev=sys.argv[4]=='1'; smin=int(sys.argv[5]); smax=int(sys.argv[6])
        S=B|fronts(R,smin,smax,betas,rev)
    Z0=np.zeros(4**W,bool); Z0[list(S)]=True
    print('R',R,mode,sys.argv[3:],'seed',len(S),'unforced',run(R,Z0,verbose=len(sys.argv)>7))
