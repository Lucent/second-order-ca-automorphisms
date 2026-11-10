/* Additive conserved quantities h: A^r -> GF(q) of F.  Unknowns h (4^r), g (4^{r+1}); for each w in A^{r+2}:
   h(F(w)) - h(w[1..r]) - g(w[1..r+1]) + g(w[0..r]) = 0.   dim_h = nullity - 1; nontrivial = dim_h - 4^{r-1}.
   Rank over GF(q), q = 32749 (nullity over GF(q) >= nullity over Q, so nontrivial=0 over GF(q) => none over Q). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define Q 32749
static int rule;
static inline int f(int a,int b,int c){return (rule>>((a<<2)|(b<<1)|c))&1;}
static inline int Floc(int s0,int s1,int s2){int x=s1>>1,y=s1&1;return (y<<1)|(x^f(s0&1,y,s2&1));}
static int64_t pw(int64_t a,int64_t e){int64_t r=1;a%=Q;while(e){if(e&1)r=r*a%Q;a=a*a%Q;e>>=1;}return r;}
int main(int argc,char**argv){ rule=atoi(argv[1]); int r=atoi(argv[2]);
  int nh=1<<(2*r), ng=1<<(2*(r+1)), nc=nh+ng, ne=1<<(2*(r+2));
  int32_t *A=calloc((size_t)ne*nc,4);
  for(int w=0;w<ne;w++){ int b[16]; int m=r+2; for(int i=0;i<m;i++) b[i]=(w>>(2*(m-1-i)))&3;
    int fw=0; for(int i=0;i<r;i++) fw=(fw<<2)|Floc(b[i],b[i+1],b[i+2]);
    int mid=0; for(int i=1;i<=r;i++) mid=(mid<<2)|b[i];
    int gr=w&(ng-1), gl=w>>2;
    int32_t *row=A+(size_t)w*nc;
    row[fw]=(row[fw]+1)%Q; row[mid]=(row[mid]+Q-1)%Q; row[nh+gr]=(row[nh+gr]+Q-1)%Q; row[nh+gl]=(row[nh+gl]+1)%Q; }
  int rank=0;
  for(int c=0;c<nc && rank<ne;c++){ int piv=-1; for(int i=rank;i<ne;i++) if(A[(size_t)i*nc+c]){piv=i;break;}
    if(piv<0) continue;
    if(piv!=rank){ for(int j=0;j<nc;j++){int32_t t=A[(size_t)piv*nc+j];A[(size_t)piv*nc+j]=A[(size_t)rank*nc+j];A[(size_t)rank*nc+j]=t;} }
    int32_t *pr=A+(size_t)rank*nc; int64_t inv=pw(pr[c],Q-2); for(int j=c;j<nc;j++) pr[j]=pr[j]*inv%Q;
    for(int i=rank+1;i<ne;i++){ int32_t *ri=A+(size_t)i*nc; int64_t fct=ri[c]; if(!fct) continue; for(int j=c;j<nc;j++) if(pr[j]) ri[j]=(int32_t)((ri[j]+(int64_t)(Q-fct)*pr[j])%Q); }
    rank++; }
  int null=nc-rank, dimh=null-1, triv=1<<(2*(r-1));
  printf("rule %d r %d: rank %d nullity %d dim_h %d trivial %d NONTRIVIAL %d%s\n",rule,r,rank,null,dimh,triv,dimh-triv,(dimh-triv)>0?"  <--":"");
  return 0;}
