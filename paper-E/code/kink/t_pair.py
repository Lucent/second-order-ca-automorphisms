import sim, numpy as np
T=80; L=2*T+20
def th(a,m): return [m if k<a else 1-m for k in range(L)]
for m in (0,1):
  for a in (6,7):
    for b in (6,7):
      PI=np.array(th(a,m),np.uint8); XI=np.array(th(b,m),np.uint8)
      cur,off,f=sim.run(94,[PI],[XI],T,(m,m))
      g,o=sim.outputs(cur,off,T,30,30)
      k=lambda s:[x for x in range(1,len(s)) if s[x]!=s[x-1]]
      print(m,a,b,'fires',sorted(f[0]),'Gamma',g[0][:3],k(g[0]),'O',o[0][:3],k(o[0]))
