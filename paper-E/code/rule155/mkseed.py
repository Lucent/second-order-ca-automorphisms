# usage: mkseed.py R betas rev smin smax genuine out.bin   (betas '-' = none)
import sys
from seeds import base, fronts, save
R=int(sys.argv[1]); b=sys.argv[2]; rev=sys.argv[3]=='1'; smin=int(sys.argv[4]); smax=int(sys.argv[5]); g=sys.argv[6]=='1'
S=base(R)
if b!='-': S|=fronts(R,smin,smax,tuple(int(c) for c in b),rev,genuine_only=g)
save(S,R,sys.argv[7]); print('seed',len(S))
