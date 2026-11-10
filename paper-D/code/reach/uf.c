#include <stdint.h>
/* union-find: on return p[i] is the smallest element of the class of i */
static int64_t find(int64_t*p,int64_t a){ while(p[a]!=a){ p[a]=p[p[a]]; a=p[a]; } return a; }
void uf(int64_t n,int64_t m,const int64_t*u,const int64_t*v,int64_t*p){
  for(int64_t i=0;i<m;i++){ int64_t a=find(p,u[i]),b=find(p,v[i]); if(a<b) p[b]=a; else if(b<a) p[a]=b; }
  for(int64_t i=0;i<n;i++) p[i]=find(p,i);
}
