/* moves.h -- exact second-order dynamics on words and the edge-move tables.
   Bit d of a word is column r-d; cells to the right of r are (0,0).
   RULE is 110, 42 or 46; T is the number of cells in a word, MSK = 2^T - 1. */
#include <stdint.h>
typedef uint64_t u64;
static int RULE, K, L, T; static u64 MSK;
static u64 fword(u64 y){ u64 a=y>>1, b=y, c=(y<<1)&MSK; u64 r;
  if(RULE==110) r=(b^c)|(b&~a);
  else if(RULE==42) r=c&~(a&b);
  else r=(c&~(a&b))|(b&~a);
  return r&MSK; }
static void fwd(u64*X,u64*Y){ u64 n=*X^fword(*Y); *X=*Y; *Y=n; }
static void bwd(u64*X,u64*Y){ u64 o=*Y^fword(*X); *Y=*X; *X=o; }
static void go(u64*X,u64*Y,int dt){ while(dt>0){fwd(X,Y);dt--;} while(dt<0){bwd(X,Y);dt++;} }
static void table(int xr,int yr,int xl,int yl,int bar,int*nx,int*ny){
  *nx=xr; *ny=yr;
  if(!bar){ if(xr==0&&yr==1){*nx=1;*ny=0;}
            else if(xr==1){ if(RULE==110) *ny=yr^xl; else if(RULE==42) *ny=yr^1^xl; } }
  else    { if(xr==1&&yr==0){*nx=0;*ny=1;}
            else if(yr==1){ if(RULE==110) *nx=xr^yl; else if(RULE==42) *nx=xr^1^yl; } }
}
/* apply the move to the time-0 state (X,Y); returns 1 if it acts */
static int move(u64*X,u64*Y,int m,int bar){
  u64 a=*X,b=*Y; go(&a,&b,m);
  int nx,ny; table(a&1,b&1,(a>>1)&1,(b>>1)&1,bar,&nx,&ny);
  if(nx==(int)(a&1)&&ny==(int)(b&1)) return 0;
  a=(a&~1ULL)|nx; b=(b&~1ULL)|ny; go(&a,&b,-m); *X=a; *Y=b; return 1; }
