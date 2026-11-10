// Exhaustive checks of the radius-7 block map h = h_* given by the tables TX, TY (files tx7.bin, ty7.bin, 2^15 bytes
// each), written at cell level
// from the definition of F = F_95 and of the table map:
//   y'_i = TY[ sum_m b_m << (m+7) ],  b_m = y_{i+m} (m odd), x_{i+m} (m even), |m|<=7
//   x'_i = TX[ sum_m b_m << (m+7) ],  b_m = x_{i+m} (m odd), y_{i+m} (m even)
// (A) h F = F h at cell 0: each output bit of both sides is a function of the 17 bits of one class in cells -8..8
//     (class 0 = {y_k, k even; x_k, k odd} for the y-bit, class 1 for the x-bit); all 2^17 values, the other bits random.
// (B) (sigma h)^2 = id at cell 0: each output bit is a function of the 29 bits of its class in cells -12..16; all 2^29.
// (C) both tests on 2^NC windows with all bits random.
// usage: verify7 tx7.bin ty7.bin [FLIP [NB [NC]]]
//   FLIP >= 0: index of a TY entry to invert before the checks (negative control); FLIP < 0: table unchanged.
//   NB (default 29): (B) enumerates only the 2^NB windows whose bits NB..28 are 0.   NC (default 26).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#define OFF 24
#define NC 49
static uint8_t TX[1<<15], TY[1<<15];
typedef struct { uint8_t x[NC], y[NC]; } cfg;      // cell k stored at index k+OFF
static inline int f95(int a,int b,int c){ return (95 >> (4*a+2*b+c)) & 1; }
static inline int gy(const cfg*c,int i){ int w=0; for(int m=-7;m<=7;m++){ int b=(m&1)? c->y[i+m+OFF] : c->x[i+m+OFF]; w|=b<<(m+7);} return TY[w]; }
static inline int gx(const cfg*c,int i){ int w=0; for(int m=-7;m<=7;m++){ int b=(m&1)? c->x[i+m+OFF] : c->y[i+m+OFF]; w|=b<<(m+7);} return TX[w]; }
static inline uint64_t rng(uint64_t*s){ uint64_t z=(*s+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31); }
static void fill_random(cfg*c,uint64_t*s){ uint64_t a=rng(s),b=rng(s); for(int k=0;k<NC;k++){ c->x[k]=(a>>k)&1; c->y[k]=(b>>k)&1; } }
// set the bits of class cl in cells lo..hi from mask
static void set_class(cfg*c,int cl,int lo,int hi,uint64_t mask){ int j=0; for(int k=lo;k<=hi;k++,j++){ int b=(mask>>j)&1; int even=((k&1)==0); if(cl==0){ if(even) c->y[k+OFF]=b; else c->x[k+OFF]=b; } else { if(even) c->x[k+OFF]=b; else c->y[k+OFF]=b; } } }
// commutation at cell 0: returns bit0 = x-bit mismatch, bit1 = y-bit mismatch
static int comm(const cfg*c){
  cfg Fc, gc;
  for(int k=-7;k<=7;k++){ Fc.x[k+OFF]=c->y[k+OFF]; Fc.y[k+OFF]=c->x[k+OFF]^f95(c->y[k-1+OFF],c->y[k+OFF],c->y[k+1+OFF]); }
  int gFx=gx(&Fc,0), gFy=gy(&Fc,0);
  for(int k=-1;k<=1;k++){ gc.x[k+OFF]=gx(c,k); gc.y[k+OFF]=gy(c,k); }
  int Fgx=gc.y[OFF], Fgy=gc.x[OFF]^f95(gc.y[-1+OFF],gc.y[OFF],gc.y[1+OFF]);
  return (gFx!=Fgx) | ((gFy!=Fgy)<<1);
}
// (sigma h)^2 at cell 0 vs identity
static int invol(const cfg*c){
  cfg s;   // s = sigma h (c): s_i = h(c)_{i+1}, needed for i in -6..8 (h reads cells 1-7..1+7 of s)
  for(int i=-6;i<=8;i++){ s.x[i+OFF]=gx(c,i+1); s.y[i+OFF]=gy(c,i+1); }
  int x2=gx(&s,1), y2=gy(&s,1);
  return (x2!=c->x[OFF]) | ((y2!=c->y[OFF])<<1);
}
int main(int argc,char**argv){
  if(argc<3){ fprintf(stderr,"usage: %s tx7.bin ty7.bin [FLIP [NB [NC]]]\n",argv[0]); return 2; }
  FILE*f=fopen(argv[1],"rb"); if(!f||fread(TX,1,1<<15,f)!=(1<<15)){ fprintf(stderr,"cannot read %s\n",argv[1]); return 1; } fclose(f);
  f=fopen(argv[2],"rb"); if(!f||fread(TY,1,1<<15,f)!=(1<<15)){ fprintf(stderr,"cannot read %s\n",argv[2]); return 1; } fclose(f);
  int nbB = 29, nbC = 26; long flip=-1;
  if(argc>3) flip=atol(argv[3]);
  if(argc>4) nbB=atoi(argv[4]);
  if(argc>5) nbC=atoi(argv[5]);
  if(flip>=(1L<<15) || nbB<0 || nbB>29 || nbC<0 || nbC>40){ fprintf(stderr,"bad arguments\n"); return 2; }
  if(flip>=0){ TY[flip]^=1; printf("NEGATIVE CONTROL: TY[%ld] flipped\n",flip); }
  // (A)
  for(int cl=0;cl<2;cl++){
    long bad=0;
    #pragma omp parallel for reduction(+:bad)
    for(long m=0;m<(1L<<17);m++){ uint64_t s=0x1234567ULL*(m+1)+cl; cfg c; fill_random(&c,&s); set_class(&c,cl,-8,8,m); int r=comm(&c); bad += (cl==0)? ((r>>1)&1) : (r&1); }
    printf("(A) class %d (%s-bit): 2^17 windows, mismatches %ld\n",cl,cl==0?"y":"x",bad);
  }
  // (B)
  for(int cl=0;cl<2;cl++){
    long bad=0; long N=1L<<nbB;
    #pragma omp parallel for reduction(+:bad) schedule(static)
    for(long m=0;m<N;m++){ uint64_t s=0x9876543ULL*(m+1)+cl; cfg c; fill_random(&c,&s); set_class(&c,cl,-12,16,(uint64_t)m); int r=invol(&c); bad += (cl==0)? ((r>>1)&1) : (r&1); }
    printf("(B) class %d (%s-bit): 2^%d windows, violations of (sigma h)^2=id: %ld\n",cl,cl==0?"y":"x",nbB,bad);
  }
  // (C)
  { long badc=0,badi=0; long N=1L<<nbC;
    #pragma omp parallel for reduction(+:badc,badi)
    for(long m=0;m<N;m++){ uint64_t s=0xABCDEF01ULL*(m+7); cfg c; fill_random(&c,&s); badc += comm(&c)!=0; badi += invol(&c)!=0; }
    printf("(C) 2^%d fully random windows: commutation mismatches %ld, involution violations %ld\n",nbC,badc,badi);
  }
  return 0;
}
