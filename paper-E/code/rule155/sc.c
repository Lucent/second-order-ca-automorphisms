// Closure that forces only the value at the point S (same closure as sonly.py), Gauss-Seidel; monotone, so the least fixed point is the same.
// usage: sc R seedfile [maxit]   seedfile = packed little-endian bits of length 4^W (seeds.save format)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int R,W,L0; static uint64_t MASK;
static inline uint64_t F(uint64_t a,uint64_t b,uint64_t c){ // rule 155 bitwise
  uint64_t r=0;
  for(int k=0;k<8;k++) if((155>>k)&1){
    uint64_t t=((k&4)?a:~a)&((k&2)?b:~b)&((k&1)?c:~c); r|=t; }
  return r;
}
static int f1(int k){return (155>>k)&1;}
int main(int argc,char**argv){
  R=atoi(argv[1]); W=2*R+1; L0=W+4; MASK=(1ULL<<W)-1;
  int maxit=argc>3?atoi(argv[3]):100000;
  uint64_t NW=1ULL<<(2*W);
  uint8_t *Z=calloc(NW,1);
  FILE*fp=fopen(argv[2],"rb"); uint64_t nb=(NW+7)/8; if(nb<8)nb=8;
  uint8_t*buf=malloc(nb); size_t got=fread(buf,1,nb,fp); fclose(fp);
  if(got<(NW+7)/8){fprintf(stderr,"short seed file\n");return 1;}
  uint64_t seed=0; for(uint64_t w=0;w<NW;w++){Z[w]=(buf[w>>3]>>(w&7))&1; seed+=Z[w];}
  // sat[k][fm]: k=yW<<2|yC<<1|yE, fm=fN<<3|fW<<2|fC<<1|fE
  uint8_t sat[8][16];
  for(int k=0;k<8;k++)for(int fm=0;fm<16;fm++){int s=0;
    for(int g=0;g<16;g++){ if(g&fm) continue; int gN=(g>>3)&1,gW=(g>>2)&1,gC=(g>>1)&1,gE=g&1;
      int kk=k^(gW<<2|gC<<1|gE); if((gN^1)==(f1(kk)^f1(k))) s=1;}
    sat[k][fm]=s;}
  uint64_t forced=seed; printf("R %d seed %llu unforced %llu\n",R,(unsigned long long)seed,(unsigned long long)(NW-forced)); fflush(stdout);
  int it;
  for(it=0;it<maxit;it++){
    uint64_t add=0;
    #pragma omp parallel for schedule(dynamic,4096) reduction(+:add)
    for(uint64_t w=0;w<NW;w++){
      if(Z[w]) continue;
      uint64_t X=w>>W, Y=w&MASK;
      for(int e=0;e<256;e++){
        uint64_t Xe=((uint64_t)(e&3)<<(W+2))|(X<<2)|((e>>2)&3);
        uint64_t Ye=((uint64_t)((e>>4)&3)<<(W+2))|(Y<<2)|((e>>6)&3);
        uint64_t r1=Xe^F(Ye>>1,Ye,Ye<<1);
        uint64_t r2=Ye^F(r1>>1,r1,r1<<1);
        uint64_t iN=(((r1>>2)&MASK)<<W)|((r2>>2)&MASK);
        uint64_t iW=(((Ye>>3)&MASK)<<W)|((r1>>3)&MASK);
        uint64_t iC=(((Ye>>2)&MASK)<<W)|((r1>>2)&MASK);
        uint64_t iE=(((Ye>>1)&MASK)<<W)|((r1>>1)&MASK);
        int k=(int)((r1>>(R+1))&7);
        int fm=(Z[iN]<<3)|(Z[iW]<<2)|(Z[iC]<<1)|Z[iE];
        if(!sat[k][fm]){Z[w]=1;add++;break;}
      }
    }
    forced+=add;
    printf(" iter %d new %llu unforced %llu\n",it,(unsigned long long)add,(unsigned long long)(NW-forced)); fflush(stdout);
    if(!add) break;
  }
  printf("RESULT R %d seed %llu unforced %llu iters %d\n",R,(unsigned long long)seed,(unsigned long long)(NW-forced),it);
  if(argc>4){FILE*o=fopen(argv[4],"wb"); fwrite(Z,1,NW,o); fclose(o);}
  return 0;
}
