/* Hidden-bit test (W1)+(H1): nodes (w,d), w in A^k, d in {1,2,3}: partner w' = w with cell j XOR d.
   Successors: for all contexts (a,b): (F(a w b), F(a w' b)) must again differ exactly at cell j (only there).
   Greatest set of nodes closed under successors = all (W, iota) candidates satisfying W1 and H1.
   Prints the number of surviving nodes per j (0 = no hidden-bit wall of this shape). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int rule,K;
static inline int f(int a,int b,int c){return (rule>>((a<<2)|(b<<1)|c))&1;}
static inline int Floc(int s0,int s1,int s2){int x=s1>>1,y=s1&1;return (y<<1)|(x^f(s0&1,y,s2&1));}
static int64_t img(int64_t w,int a,int b){int buf[24]; buf[0]=a; for(int i=0;i<K;i++)buf[1+i]=(w>>(2*(K-1-i)))&3; buf[K+1]=b;
  int64_t e=0; for(int i=0;i<K;i++) e=(e<<2)|Floc(buf[i],buf[i+1],buf[i+2]); return e;}
int main(int argc,char**argv){ rule=atoi(argv[1]); K=atoi(argv[2]); int64_t N=1LL<<(2*K);
  int64_t *im=malloc(N*16*8); for(int64_t w=0;w<N;w++) for(int c=0;c<16;c++) im[w*16+c]=img(w,c>>2,c&3);
  int total=0;
  for(int j=1;j<=K-2;j++){ int sh=2*(K-1-j); int64_t M=N*3;
    uint8_t *bad=calloc(M,1); int32_t *succ=malloc(M*16*4); int32_t *indeg=calloc(M+1,4);
    for(int64_t n=0;n<M;n++){ int64_t w=n/3; int d=n%3+1; int64_t w2=w^((int64_t)d<<sh);
      for(int c=0;c<16;c++){ int64_t u=im[w*16+c], u2=im[w2*16+c], x=u^u2;
        if(x==0 || (x & ~(3LL<<sh))){ bad[n]=1; succ[n*16+c]=-1; } else { int64_t m=u*3+((x>>sh)&3)-1; succ[n*16+c]=m; } } }
    /* reverse adjacency (CSR) */
    for(int64_t n=0;n<M;n++) for(int c=0;c<16;c++){int32_t m=succ[n*16+c]; if(m>=0) indeg[m+1]++;}
    for(int64_t n=0;n<M;n++) indeg[n+1]+=indeg[n];
    int32_t *pred=malloc((size_t)indeg[M]*4+4), *pos=malloc(M*4); for(int64_t n=0;n<M;n++)pos[n]=indeg[n];
    for(int64_t n=0;n<M;n++) for(int c=0;c<16;c++){int32_t m=succ[n*16+c]; if(m>=0) pred[pos[m]++]=n;}
    int32_t *q=malloc(M*4); int64_t qh=0,qt=0; for(int64_t n=0;n<M;n++) if(bad[n]) q[qt++]=n;
    while(qh<qt){ int32_t m=q[qh++]; for(int32_t e=indeg[m];e<indeg[m+1];e++){int32_t p=pred[e]; if(!bad[p]){bad[p]=1;q[qt++]=p;}} }
    int64_t surv=0, ex=-1; for(int64_t n=0;n<M;n++) if(!bad[n]){surv++; if(ex<0)ex=n;}
    printf("rule %d k %d j %d: surviving (word,flip) nodes %lld",rule,K,j,(long long)surv);
    if(surv){ printf("  e.g. w="); for(int i=0;i<K;i++) printf("%d",(int)(((ex/3)>>(2*(K-1-i)))&3)); printf(" d=%d  <--",(int)(ex%3+1)); }
    printf("\n"); total+=surv>0;
    free(bad);free(succ);free(indeg);free(pred);free(pos);free(q); }
  return 0;}
