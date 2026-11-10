/* (W1) search: graph on A^k, edge w -> (sigma^s F^p c)[0,k) for every c with c[0,k)=w (all 4^{2p} contexts).
   A proper set W subset A^k with "0 in Occ(c) => 0 in Occ(sigma^s F^p c)" exists iff there is a sink SCC != A^k.
   Prints the number of sink SCCs, the size of the smallest one and up to 12 of them (size, a member) when the
   graph is not a single closed class. A word fixed by sigma^s F^p in every context is a sink SCC of size 1.
   Usage: w1 RULE K P S */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int rule,K,P,S; static int64_t N; static int NC;
static inline int f(int a,int b,int c){return (rule>>((a<<2)|(b<<1)|c))&1;}
static inline int Floc(int s0,int s1,int s2){int x=s1>>1,y=s1&1;return (y<<1)|(x^f(s0&1,y,s2&1));}
static int L,R; /* extra cells left/right */
static int64_t succ(int64_t w,int ctx){
  int buf[80]; int m=L+K+R;
  for(int i=0;i<L;i++) buf[i]=(ctx>>(2*i))&3;
  for(int i=0;i<K;i++) buf[L+i]=(w>>(2*(K-1-i)))&3;
  for(int i=0;i<R;i++) buf[L+K+i]=(ctx>>(2*(L+i)))&3;
  int len=m; for(int t=0;t<P;t++){for(int i=0;i<len-2;i++)buf[i]=Floc(buf[i],buf[i+1],buf[i+2]);len-=2;}
  /* after P steps buf[j] is cell (j - L + P) of F^P c ; want cells s..s+K-1 -> j = s+L-P */
  int off=S+L-P; int64_t e=0; for(int i=0;i<K;i++) e=(e<<2)|buf[off+i]; return e;}
int main(int argc,char**argv){
  rule=atoi(argv[1]);K=atoi(argv[2]);P=atoi(argv[3]);S=atoi(argv[4]);
  L=P-S; R=P+S; if(L<0||R<0){fprintf(stderr,"need |s|<=p\n");return 1;}
  N=1LL<<(2*K); NC=1<<(2*(L+R));
  int32_t *idx=malloc(N*4),*low=malloc(N*4),*comp=malloc(N*4),*stk=malloc(N*4),*cs=malloc(N*4),*cc=malloc(N*4);
  uint8_t *on=calloc(N,1);
  for(int64_t i=0;i<N;i++)idx[i]=-1;
  int32_t counter=0,sp=0,ncomp=0;
  for(int64_t root=0;root<N;root++){ if(idx[root]>=0)continue;
    int csp=0; cs[csp]=root; cc[csp]=0; csp++; idx[root]=low[root]=counter++; stk[sp++]=root; on[root]=1;
    while(csp){ int32_t v=cs[csp-1];
      if(cc[csp-1]<NC){ int32_t w=succ(v,cc[csp-1]++);
        if(idx[w]<0){ idx[w]=low[w]=counter++; stk[sp++]=w; on[w]=1; cs[csp]=w; cc[csp]=0; csp++; }
        else if(on[w] && idx[w]<low[v]) low[v]=idx[w];
      } else { csp--;
        if(low[v]==idx[v]){ int32_t w; do{ w=stk[--sp]; on[w]=0; comp[w]=ncomp; }while(w!=v); ncomp++; }
        if(csp){ int32_t u=cs[csp-1]; if(low[v]<low[u]) low[u]=low[v]; }
      } } }
  uint8_t *nonsink=calloc(ncomp,1); int64_t *size=calloc(ncomp,8); int64_t *rep=malloc(ncomp*8);
  for(int64_t v=0;v<N;v++){ size[comp[v]]++; rep[comp[v]]=v; for(int c=0;c<NC;c++) if(comp[succ(v,c)]!=comp[v]){nonsink[comp[v]]=1;break;} }
  int nsink=0; int64_t minsink=N; for(int c=0;c<ncomp;c++) if(!nonsink[c]){ nsink++; if(size[c]<minsink) minsink=size[c]; }
  printf("rule %d k %d p %d s %d: %d SCCs, %d sink SCC(s)",rule,K,P,S,ncomp,nsink);
  if(nsink==1 && ncomp==1) printf(" -> A^k is the only closed set\n");
  else { printf(" -> PROPER CLOSED SETS (smallest of size %lld):",(long long)minsink); int shown=0; for(int c=0;c<ncomp&&shown<12;c++) if(!nonsink[c]){ printf(" [size %lld rep ",(long long)size[c]); for(int i=0;i<K;i++) printf("%d",(int)((rep[c]>>(2*(K-1-i)))&3)); printf("]"); shown++;} printf("\n"); }
  return 0;}
