import sim, numpy as np
T=80; L=2*T+20
z=lambda s:[x for x in range(len(s)) if s[x]==0]
for a in (8,9):
  for b in (8,9):
    PI=np.ones(L,np.uint8); PI[a]=0; XI=np.ones(L,np.uint8); XI[b]=0
    cur,off,f=sim.run(95,[PI],[XI],T,(1,1)); G,O=sim.outputs(cur,off,T,40,40)
    print(a,b,'fires',sorted((int(x),int(y)) for x,y in f[0]),'Gamma zeros',z(G[0]),'O zeros',z(O[0]))
