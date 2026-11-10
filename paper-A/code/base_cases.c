/* base_cases.c -- orbits of the edge moves on a finite window.

   Rules 110, 42, 46 (second-order, F(x,y)_i = (y_i, x_i ^ f(y_{i-1},y_i,y_{i+1}))).
   A configuration c is 0-tailed with edge r: c_i = (0,0) for i > r, c_r != (0,0).
   Representation: bit d of a word is column r-d.  A state is a pair of words
   (X,Y) = (row t-1, row t); the state at time 0 is c itself, c_i = (X_d, Y_d).

   Moves (Lemmas cl:lem:edgemove, cl:lem:moves1446), at time m: run F^m, change the
   edge cell (x_r,y_r) of the state at time m by the move table, run F^{-m}.
     110:  M : 01->10;  if x_r=1: (1,y)->(1, y ^ x_{r-1})
           Mb: 10->01;  if y_r=1: (x,1)->(x ^ y_{r-1}, 1)
      42:  M : 01->10;  if x_r=1: (1,y)->(1, y ^ 1 ^ x_{r-1})
           Mb: 10->01;  if y_r=1: (x,1)->(x ^ 1 ^ y_{r-1}, 1)
      46:  M : 01->10 only;  Mb: 10->01 only.

   Mode "cone": the window is W_k = [r-3k+1, r] (L = 3k cells), every value of it with
   edge cell != (0,0) is a vertex (3*4^(3k-1) vertices), and the columns left of the
   window are padding.  By the light-cone bound (Lemma cl:lem:movelocal), whether a move
   at time |m| <= k-1 acts and its result depend only on c on [r-3|m|, r], a subset of W_k,
   and the move changes only [r-|m|, r].  The program checks this: every move is computed
   with NPAD different paddings (zeros, ones, pseudo-random) and the results on the whole
   window must coincide; cells outside [r-|m|, r] must be unchanged; the edge cell must stay
   != (0,0); every move must be undone by a move at the same time.
   The connected components (union-find) are then compared with the statement:
     110: class(w) = { w' : w' = w outside [r-k+1,r], w'_r != 00, and, if w_{r-k} = 00,
                       [w'_{r-k+1} = 11] = [w_{r-k+1} = 11] }          (F_k, Thm cl:thm:Fk)
      42: exterior != 11: class = all w' agreeing outside [r-k+1,r]       ((A) of cl:thm:Fk42)
          exterior  = 11: class = all w' with the same length of the run of 11 cells
                          starting at the tip (this equality implies (B))
      46: class(w) = { w' agreeing outside [r-k+1,r], w'_r != 00, same bond bits
                       chi_i, r-k <= i <= r-1 }                          (Thm cl:thm:Fk46)
   Mode "short" (confirmation only): window [r-k, r] and the columns left of it filled with
   zeros or with pseudo-random cells, one run per filling.

   Usage: base_cases RULE K MODE      (MODE = cone | short)
   Prints per exterior cell the number of classes, and "ALL CHECKS PASSED" or exits 1. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "moves.h"
static u64 rs=88172645463325252ULL;
static u64 rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static int *par;
static int fnd(int x){ while(par[x]!=x){ par[x]=par[par[x]]; x=par[x]; } return x; }
static void uni(int a,int b){ a=fnd(a); b=fnd(b); if(a!=b) par[a]=b; }
static int cell(u64 X,u64 Y,int d){ return (int)(((X>>d)&1)<<1 | ((Y>>d)&1)); }
static int cmp(const void*x,const void*y){ const u64*a=x,*b=y; return a[0]<b[0]?-1:a[0]>b[0]?1:(a[1]<b[1]?-1:a[1]>b[1]); }
static void fail(const char*s){ printf("FAILURE: %s\n",s); exit(1); }
/* key of the claimed class; window bits X,Y (L bits) */
static u64 key(u64 X,u64 Y){
  u64 out=(X>>K)|((Y>>K)<<L);           /* cells outside [r-k+1,r] */
  int ext=cell(X,Y,K); u64 extra=0;
  if(RULE==110){ if(ext==0) extra=(cell(X,Y,K-1)==3); }
  else if(RULE==42){ if(ext==3){ int n=0; for(int d=K-1;d>=0&&cell(X,Y,d)==3;d--) n++; extra=n; } }
  else { for(int d=K;d>=1;d--){ int chi=(cell(X,Y,d)==0)^(cell(X,Y,d-1)==3); extra=(extra<<1)|chi; } }
  return (out<<8)|extra; }
int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: base_cases RULE K cone|short\n"); return 2; }
  RULE=atoi(argv[1]); K=atoi(argv[2]); int cone=!strcmp(argv[3],"cone");
  /* rule table check against the Wolfram number */
  for(int t=0;t<8;t++){ MSK=7; T=3; if((int)((fword((u64)t)>>1)&1)!=((RULE>>t)&1)) fail("rule table"); }
  L = cone ? 3*K : K+1;
  int P = 2*K+6; T=L+P; if(T>62) fail("window too wide");
  MSK=(1ULL<<T)-1; u64 WM=(1ULL<<L)-1;
  int NPAD = cone ? 4 : 1, NFILL = cone ? 1 : 6;
  long NS=1L<<(2*L);
  par=malloc(sizeof(int)*NS);
  long totviol=0;
  for(int fill=0; fill<NFILL; fill++){
    u64 padX[8],padY[8];
    for(int p=0;p<NPAD;p++){
      u64 px,py; int mode = cone ? p : (fill==0?0:2);
      if(mode==0){px=0;py=0;} else if(mode==1){px=MSK;py=MSK;} else {px=rnd();py=rnd();}
      padX[p]=px&MSK&~WM; padY[p]=py&MSK&~WM; }
    for(long s=0;s<NS;s++) par[s]=(int)s;
    long nvert=0, nedge=0;
    for(long s=0;s<NS;s++){
      u64 X=s&WM, Y=(u64)s>>L; if(!((X|Y)&1)) continue; nvert++;
      for(int m=-(K-1);m<=K-1;m++) for(int bar=0;bar<2;bar++){
        long tgt=-1; int acts=-1;
        for(int p=0;p<NPAD;p++){
          u64 a=X|padX[p], b=Y|padY[p]; int ac=move(&a,&b,m,bar);
          u64 ch=((a^(X|padX[p]))|(b^(Y|padY[p])));
          int am=m<0?-m:m;
          if((ch & ~((1ULL<<(am+1))-1))) fail("move changed a cell outside [r-|m|,r]");
          if(!((a|b)&1)) fail("edge cell became (0,0)");
          long t2=(long)((a&WM)|((b&WM)<<L));
          if(p==0){tgt=t2;acts=ac;} else if(t2!=tgt||ac!=acts) fail("result depends on the padding");
        }
        if(!acts) continue;
        /* reversibility: some move at time m maps the image back */
        { int back=0; for(int bb=0;bb<2;bb++){ u64 a=(tgt&WM)|padX[0], b=((u64)tgt>>L)|padY[0];
            move(&a,&b,m,bb); if((long)((a&WM)|((b&WM)<<L))==s) back=1; }
          if(!back) fail("move not undone by a move at the same time"); }
        nedge++; uni((int)s,(int)tgt);
      }
    }
    /* compare components with the claimed partition: bijection root <-> key */
    /* hash maps via sorting pairs */
    u64 *pr=malloc(sizeof(u64)*2*nvert); long n=0;
    for(long s=0;s<NS;s++){ u64 X=s&WM,Y=(u64)s>>L; if(!((X|Y)&1)) continue;
      pr[2*n]=(u64)fnd((int)s); pr[2*n+1]=key(X,Y); n++; }
    /* check: each root has one key, each key one root */
    long nroots=0,nkeys=0,bad=0;
    { /* root -> key */
      u64 *rk=malloc(sizeof(u64)*NS); char*seen=calloc(NS,1);
      for(long i=0;i<n;i++){ u64 r=pr[2*i]; if(!seen[r]){seen[r]=1;rk[r]=pr[2*i+1];nroots++;} else if(rk[r]!=pr[2*i+1]) bad++; }
      free(rk); free(seen); }
    { /* key -> root: sort keys */
      u64 *kk=malloc(sizeof(u64)*2*n);
      for(long i=0;i<n;i++){ kk[2*i]=pr[2*i+1]; kk[2*i+1]=pr[2*i]; }
      qsort(kk,n,2*sizeof(u64),cmp);
      for(long i=0;i<n;i++){ if(i==0||kk[2*i]!=kk[2*i-2]) nkeys++; else if(kk[2*i+1]!=kk[2*i-1]) bad++; }
      free(kk); }
    /* classes per exterior cell, per fixed cells left of the exterior */
    long cnt[4]={0,0,0,0};
    { char*seen=calloc(NS,1);
      for(long s=0;s<NS;s++){ u64 X=s&WM,Y=(u64)s>>L; if(!((X|Y)&1)) continue;
        if(((X>>(K+1))|(Y>>(K+1)))!=0) continue;           /* cells left of the exterior = (0,0) */
        int r=fnd((int)s); if(!seen[r]){ seen[r]=1; cnt[cell(X,Y,K)]++; } }
      free(seen); }
    printf("rule %d k=%d mode=%s fill=%d: vertices %ld, edges %ld, classes %ld, claimed classes %ld, mismatches %ld;"
           " classes per exterior cell 00,01,10,11 (cells further left (0,0)): %ld %ld %ld %ld\n",
           RULE,K,cone?"cone":"short",fill,nvert,nedge,nroots,nkeys,bad,cnt[0],cnt[1],cnt[2],cnt[3]);
    totviol+=bad+(nroots!=nkeys);
    free(pr);
  }
  if(totviol) fail("partition differs from the statement");
  printf("ALL CHECKS PASSED\n"); return 0; }
