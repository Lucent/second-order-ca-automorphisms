/* Window closure for F_f, f any elementary rule (the closure of ../closure/dc.c, with OpenMP and a base read from a file).
 * Window index: idx=(Xw<<W)|Yw, W=2R+1, cell k (k=0 leftmost, position k-R) at bit W-1-k.
 * Diamond: the spacetime points (t,i) with |t|+|i|<=D; the target is (0,0).  An extension of the window by D cells per
 * side at time 0 determines the windows at all points of the diamond.
 * usage: lc dc RULE R 2 BASEFILE
 * Mode dc, D=2: a window is forced if, for some extension and one of the five instances of the spacetime relation
 * through the centre, the relation has no solution with difference 1 at the centre and 0 at those of its other four
 * points whose windows are already forced.  Each iteration reads the set as at the start of the iteration.
 * The other modes of the program (full, five, sets, sound, minsup, dfs) are not used for the paper. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <omp.h>
typedef uint64_t u64;
static int RULE,R,W,D,L0,NP,NCTX; static u64 mW,NW;
static int PT[64],PI[64]; static int PIDX[20][20]; /* point index by (t+D,i+D); centre = -1 */
static u64 **MS; static int *NMS; /* minimal supports per ctx */
static inline int fb(int a,int b,int c){ return (RULE>>((a<<2)|(b<<1)|c))&1; }
static inline u64 fpl(u64 a,u64 b,u64 c){ u64 r=0; for(int k=0;k<8;k++) if((RULE>>k)&1)
  r|=((k&4)?a:~a)&((k&2)?b:~b)&((k&1)?c:~c); return r; }
static inline void fw(u64 X,u64 Y,int n,u64*Xo,u64*Yo){ u64 m=(1ULL<<(n-2))-1; u64 fy=fpl(Y>>1,Y,Y<<1);
  *Xo=(Y>>1)&m; *Yo=((X^fy)>>1)&m; }
static inline void bw(u64 X,u64 Y,int n,u64*Xo,u64*Yo){ u64 m=(1ULL<<(n-2))-1; u64 fx=fpl(X>>1,X,X<<1);
  *Xo=((Y^fx)>>1)&m; *Yo=(X>>1)&m; }
static inline int getb(const u64*Z,u64 i){ return (Z[i>>6]>>(i&63))&1; }
static double now(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static size_t nwords(){ return NW<64?1:NW/64; }
static u64 popc(const u64*Z){ u64 s=0; for(size_t i=0;i<nwords();i++) s+=__builtin_popcountll(Z[i]); return s; }
/* ---------- diamond geometry ---------- */
static void geom(){ NP=0; for(int t=-D;t<=D;t++) for(int i=-D;i<=D;i++){ int a=abs(t)+abs(i); PIDX[t+D][i+D]=-2;
  if(a<=D){ if(t==0&&i==0){PIDX[D][D]=-1;continue;} PT[NP]=t;PI[NP]=i;PIDX[t+D][i+D]=NP++; } } NCTX=1<<(4*D); }
/* ctx bits: Y(0,i) i=-D..D -> bits 0..2D ; Y(-1,i) i=-(D-1)..D-1 -> bits 2D+1.. */
static void diamondY(int ctx,int Y[20][20]){ /* Y[t+D][i+D] */
  for(int i=-D;i<=D;i++) Y[D][i+D]=(ctx>>(i+D))&1;
  for(int i=-(D-1);i<=D-1;i++) Y[D-1][i+D]=(ctx>>(2*D+1+i+D-1))&1;
  for(int t=0;t<D;t++) for(int i=-(D-1-t);i<=D-1-t;i++)
    Y[t+1+D][i+D]=Y[t-1+D][i+D]^fb(Y[t+D][i-1+D],Y[t+D][i+D],Y[t+D][i+1+D]);
  for(int t=-1;t>-D;t--) for(int i=-(D-1+t);i<=D-1+t;i++)   /* |t|=-t: stencil at (t,i), |i|<=D-1-|t| */
    Y[t-1+D][i+D]=Y[t+1+D][i+D]^fb(Y[t+D][i-1+D],Y[t+D][i+D],Y[t+D][i+1+D]);
}
static void minsup(){
  MS=malloc(sizeof(u64*)*NCTX); NMS=calloc(NCTX,sizeof(int));
  int nfree=4*D; 
  #pragma omp parallel for schedule(dynamic,16)
  for(int ctx=0;ctx<NCTX;ctx++){ int Y[20][20],G[20][20],Yp[20][20]; diamondY(ctx,Y);
    u64 *sup=malloc(sizeof(u64)<<(nfree-1)); int ns=0;
    for(int gf=0;gf<(1<<nfree);gf++){ if(!((gf>>D)&1)) continue; /* G(0,0)=bit D */
      memset(G,0,sizeof G);
      for(int i=-D;i<=D;i++) G[D][i+D]=(gf>>(i+D))&1;
      for(int i=-(D-1);i<=D-1;i++) G[D-1][i+D]=(gf>>(2*D+1+i+D-1))&1;
      for(int t=-1;t<=0;t++) for(int i=-D;i<=D;i++) Yp[t+D][i+D]=Y[t+D][i+D]^G[t+D][i+D];
      for(int t=0;t<D;t++){ for(int i=-(D-1-t);i<=D-1-t;i++){
          int d=fb(Yp[t+D][i-1+D],Yp[t+D][i+D],Yp[t+D][i+1+D])^fb(Y[t+D][i-1+D],Y[t+D][i+D],Y[t+D][i+1+D]);
          G[t+1+D][i+D]=G[t-1+D][i+D]^d; Yp[t+1+D][i+D]=Y[t+1+D][i+D]^G[t+1+D][i+D]; } }
      for(int t=-1;t>-D;t--){ for(int i=-(D-1+t);i<=D-1+t;i++){
          int d=fb(Yp[t+D][i-1+D],Yp[t+D][i+D],Yp[t+D][i+1+D])^fb(Y[t+D][i-1+D],Y[t+D][i+D],Y[t+D][i+1+D]);
          G[t-1+D][i+D]=G[t+1+D][i+D]^d; Yp[t-1+D][i+D]=Y[t-1+D][i+D]^G[t-1+D][i+D]; } }
      u64 s=0; for(int p=0;p<NP;p++) if(G[PT[p]+D][PI[p]+D]) s|=1ULL<<p;
      sup[ns++]=s; }
    /* minimal elements */
    for(int a=0;a<ns;a++) for(int b=a+1;b<ns;b++) if(__builtin_popcountll(sup[b])<__builtin_popcountll(sup[a])){u64 t=sup[a];sup[a]=sup[b];sup[b]=t;}
    u64 *m=malloc(8*ns); int nm=0;
    for(int a=0;a<ns;a++){ int dom=0; for(int b=0;b<nm;b++) if((m[b]&sup[a])==m[b]){dom=1;break;} if(!dom) m[nm++]=sup[a]; }
    MS[ctx]=realloc(m,8*(nm?nm:1)); NMS[ctx]=nm; free(sup); }
}
static inline int implied(int ctx,u64 zm){ const u64*m=MS[ctx]; int n=NMS[ctx];
  for(int k=0;k<n;k++) if(!(m[k]&zm)) return 0; return 1; }
/* ---------- patches ---------- */
/* given extended row-0 planes (length L0), compute the window index at every diamond point and the ctx */
static inline void patch(u64 X0,u64 Y0,u64*widx,int*ctx){
  u64 RX[20],RY[20]; RX[D]=X0; RY[D]=Y0;
  for(int t=1;t<=D;t++){ fw(RX[D+t-1],RY[D+t-1],L0-2*(t-1),&RX[D+t],&RY[D+t]); bw(RX[D-t+1],RY[D-t+1],L0-2*(t-1),&RX[D-t],&RY[D-t]); }
  for(int p=-1;p<NP;p++){ int t= p<0?0:PT[p], i= p<0?0:PI[p]; int Lt=L0-2*abs(t); int a=i+D-abs(t); int sh=Lt-W-a;
    widx[p+1]=(((RX[D+t]>>sh)&mW)<<W)|((RY[D+t]>>sh)&mW); }
  int c=0; for(int i=-D;i<=D;i++) c|=((Y0>>(L0-1-(i+R+D)))&1)<<(i+D);
  for(int i=-(D-1);i<=D-1;i++) c|=((X0>>(L0-1-(i+R+D)))&1)<<(2*D+1+i+D-1);
  *ctx=c; }
/* five original rules on the patch: point bits of zm; Y from ctx */
static int PB(int t,int i){ return PIDX[t+D][i+D]; }
static inline int zp(u64 zm,int t,int i){ return (zm>>PB(t,i))&1; }
static inline int five(int ctx,u64 zm){ int Y[20][20]; diamondY(ctx,Y);
  if(zp(zm,-2,0)&&zp(zm,-1,-1)&&zp(zm,-1,0)&&zp(zm,-1,1)) return 1;
  if(zp(zm,2,0)&&zp(zm,1,-1)&&zp(zm,1,0)&&zp(zm,1,1)) return 1;
  { int b=Y[D][1+D],c=Y[D][2+D]; if(fb(0,b,c)!=fb(1,b,c) && zp(zm,1,1)&&zp(zm,-1,1)&&zp(zm,0,1)&&zp(zm,0,2)) return 1; }
  { int a=Y[D][-1+D],c=Y[D][1+D]; if(fb(a,0,c)!=fb(a,1,c) && zp(zm,1,0)&&zp(zm,-1,0)&&zp(zm,0,-1)&&zp(zm,0,1)) return 1; }
  { int a=Y[D][-2+D],b=Y[D][-1+D]; if(fb(a,b,0)!=fb(a,b,1) && zp(zm,1,-1)&&zp(zm,-1,-1)&&zp(zm,0,-2)&&zp(zm,0,-1)) return 1; }
  return 0; }
/* single-stencil don't-care rules: some stencil through the target, all of whose other non-forced points are free,
   forces G(0,0)=0 by the one 5-point relation alone */
static inline int dc(int ctx,u64 zm){ int Y[20][20]; diamondY(ctx,Y);
  static const int C[5][2]={{0,0},{1,0},{-1,0},{0,1},{0,-1}};
  for(int s=0;s<5;s++){ int st=C[s][0],si=C[s][1]; int pt[5][2]={{st+1,si},{st-1,si},{st,si-1},{st,si},{st,si+1}};
    int unk[5],nu=0,q=-1;
    for(int k=0;k<5;k++){ if(pt[k][0]==0&&pt[k][1]==0){q=k;continue;} if(!zp(zm,pt[k][0],pt[k][1])) unk[nu++]=k; }
    int ok=1;
    for(int a=0;a<(1<<nu)&&ok;a++){ int g[5]={0,0,0,0,0}; g[q]=1; for(int u=0;u<nu;u++) g[unk[u]]=(a>>u)&1;
      int y[5]; for(int k=0;k<5;k++) y[k]=Y[pt[k][0]+D][pt[k][1]+D];
      int lhs=g[0]^g[1]; int rhs=fb(y[2]^g[2],y[3]^g[3],y[4]^g[4])^fb(y[2],y[3],y[4]);
      if(lhs==rhs) ok=0; }
    if(ok) return 1; }
  return 0; }
/* direct search: is there G on the diamond with G(0,0)=1, G=0 on zm, satisfying all stencils?  (same criterion as implied()) */
static int EQ=0;
static inline int dfs(int ctx,u64 zm,const u64*wid){ int Y[20][20]; diamondY(ctx,Y);
  int fr[40][2],nf=0; 
  for(int i=-D;i<=D;i++) if(i!=0 && !zp(zm,0,i)){fr[nf][0]=0;fr[nf][1]=i;nf++;}
  for(int i=-(D-1);i<=D-1;i++) if(!zp(zm,-1,i)){fr[nf][0]=-1;fr[nf][1]=i;nf++;}
  int G[20][20],Yp[20][20];
  for(long a=0;a<(1L<<nf);a++){ memset(G,0,sizeof G); G[D][D]=1; for(int k=0;k<nf;k++) G[fr[k][0]+D][fr[k][1]+D]=(a>>k)&1;
    for(int t=-1;t<=0;t++) for(int i=-D;i<=D;i++) Yp[t+D][i+D]=Y[t+D][i+D]^G[t+D][i+D];
    int bad=0;
    for(int t=0;t<D&&!bad;t++) for(int i=-(D-1-t);i<=D-1-t;i++){
      int d=fb(Yp[t+D][i-1+D],Yp[t+D][i+D],Yp[t+D][i+1+D])^fb(Y[t+D][i-1+D],Y[t+D][i+D],Y[t+D][i+1+D]);
      int g=G[t-1+D][i+D]^d; if(g&&zp(zm,t+1,i)){bad=1;break;} G[t+1+D][i+D]=g; Yp[t+1+D][i+D]=Y[t+1+D][i+D]^g; }
    for(int t=-1;t>-D&&!bad;t--) for(int i=-(D-1+t);i<=D-1+t;i++){
      int d=fb(Yp[t+D][i-1+D],Yp[t+D][i+D],Yp[t+D][i+1+D])^fb(Y[t+D][i-1+D],Y[t+D][i+D],Y[t+D][i+1+D]);
      int g=G[t+1+D][i+D]^d; if(g&&zp(zm,t-1,i)){bad=1;break;} G[t-1+D][i+D]=g; Yp[t-1+D][i+D]=Y[t-1+D][i+D]^g; }
    if(!bad && EQ){ /* points with identical windows must carry identical G */
      for(int p=-1;p<NP&&!bad;p++) for(int q=p+1;q<NP;q++) if(wid[p+1]==wid[q+1]){
        int tp=p<0?0:PT[p], ip=p<0?0:PI[p];
        if(G[tp+D][ip+D]!=G[PT[q]+D][PI[q]+D]){bad=1;break;} } }
    if(!bad) return 0; }
  return 1; }
/* restricted rule sets */
static int NRS; static u64 *RS; static u64 **RC; /* RC[k] = bitmap over ctx */
static inline int sets(int ctx,u64 zm){ for(int k=0;k<NRS;k++) if((RS[k]&zm)==RS[k] && ((RC[k][ctx>>6]>>(ctx&63))&1)) return 1; return 0; }
static void loadrules(const char*fn){ FILE*fp=fopen(fn,"r"); if(!fp){perror(fn);exit(2);} NRS=0; RS=malloc(8*100000);
  unsigned long long v; while(fscanf(fp,"%llx",&v)==1) RS[NRS++]=v; fclose(fp);
  RC=malloc(sizeof(u64*)*NRS);
  for(int k=0;k<NRS;k++){ RC[k]=calloc(NCTX/64+1,8); for(int c=0;c<NCTX;c++) if(implied(c,RS[k])) RC[k][c>>6]|=1ULL<<(c&63); } }
static int MODE; /* 0 full 1 five 2 sets 3 dc */
static unsigned char *DCT=NULL;
static inline int crit(int ctx,u64 zm,const u64*wid){ if(DCT){ u64 k=((u64)ctx<<NP)|zm; return (DCT[k>>3]>>(k&7))&1; } return MODE==0?implied(ctx,zm):(MODE==1?five(ctx,zm):(MODE==2?sets(ctx,zm):(MODE==3?dc(ctx,zm):dfs(ctx,zm,wid)))); }
static void setup(int rule,int r,int d){ RULE=rule;R=r;D=d;W=2*R+1;L0=W+2*D;mW=(1ULL<<W)-1;NW=1ULL<<(2*W); geom(); }
static inline void ext(u64 w,u64 e,u64*X0,u64*Y0){ u64 Xc=w>>W,Yc=w&mW; u64 ml=(1ULL<<D)-1;
  u64 xl=e&ml, xr=(e>>D)&ml, yl=(e>>(2*D))&ml, yr=(e>>(3*D))&ml;
  *X0=(xl<<(W+D))|(Xc<<D)|xr; *Y0=(yl<<(W+D))|(Yc<<D)|yr; }
static int closure(const char*basefn,const char*stampfn){
  size_t nw=nwords(); u64*Z=calloc(nw,8),*Zn=calloc(nw,8); uint16_t*ST=calloc(NW,2); uint32_t*WX=calloc(NW,4);
  FILE*fp=fopen(basefn,"rb"); if(!fp){perror(basefn);return 2;} size_t got=fread(Z,8,nw,fp); fclose(fp); (void)got;
  if(getenv("LC_INIT")){ fp=fopen(getenv("LC_INIT"),"rb"); if(!fp){perror("init");return 2;} u64*T=calloc(nw,8); got=fread(T,8,nw,fp); fclose(fp);
    for(size_t i=0;i<nw;i++) Z[i]|=T[i]; free(T); }
  if(NW<64) Z[0]&=(1ULL<<NW)-1;
  u64*AL=NULL; if(getenv("LC_ALLOW")){ fp=fopen(getenv("LC_ALLOW"),"rb"); AL=calloc(nw,8); got=fread(AL,8,nw,fp); fclose(fp); }
  int gs=getenv("LC_GS")!=NULL; u64 cnt=popc(Z); u64 NE=1ULL<<(4*D);
  printf("rule %d R %d D %d mode %d windows %llu base %llu points %d\n",RULE,R,D,MODE,(unsigned long long)NW,(unsigned long long)cnt,NP); fflush(stdout);
  double t0=now(); int it=0;
  while(1){ it++; u64 before=cnt; memcpy(Zn,Z,nw*8); u64*Zr=gs?Zn:Z;
    #pragma omp parallel for schedule(dynamic,4)
    for(size_t wi=0;wi<nw;wi++){ u64 word=Z[wi]; if(NW>=64&&word==~0ULL) continue;
      for(int b=0;b<64&&(u64)(wi*64+b)<NW;b++){ if((word>>b)&1) continue; u64 w=wi*64+b; if(AL&&!((AL[wi]>>b)&1)) continue;
        for(u64 e=0;e<NE;e++){ u64 X0,Y0,wid[64]; int ctx; ext(w,e,&X0,&Y0); patch(X0,Y0,wid,&ctx);
          u64 zm=0; for(int p=0;p<NP;p++) if(getb(Zr,wid[p+1])) zm|=1ULL<<p;
          if(crit(ctx,zm,wid)){ __atomic_fetch_or(&Zn[wi],1ULL<<b,__ATOMIC_RELAXED); ST[w]=it; WX[w]=(uint32_t)e; break; } } } }
    u64*t=Z;Z=Zn;Zn=t; cnt=popc(Z);
    printf(" iter %d %llu (%.1f s)\n",it,(unsigned long long)cnt,now()-t0); fflush(stdout);
    if(cnt==before) break; }
  printf("RESULT rule %d R %d D %d mode %d unforced %llu iterations %d time %.1f\n",RULE,R,D,MODE,(unsigned long long)(NW-cnt),it-1,now()-t0);
  if(stampfn){ char fn[600]; snprintf(fn,600,"%s.stamp",stampfn); fp=fopen(fn,"wb"); fwrite(ST,2,NW,fp); fclose(fp);
    snprintf(fn,600,"%s.ext",stampfn); fp=fopen(fn,"wb"); fwrite(WX,4,NW,fp); fclose(fp);
    snprintf(fn,600,"%s.Z",stampfn); fp=fopen(fn,"wb"); fwrite(Z,8,nw,fp); fclose(fp); }
  return NW==cnt?0:1; }
/* ---------- soundness ---------- */
static u64*G;
static void gmap(int a,int b,int cst){ size_t nw=nwords(); memset(G,0,nw*8);
  for(u64 w=0;w<NW;w++){ u64 X=w>>W,Y=w&mW; int n=W; int yh;
    if(cst>=0) yh=cst&1; else { for(int s=0;s<abs(b);s++){ u64 X1,Y1; if(b>0) fw(X,Y,n,&X1,&Y1); else bw(X,Y,n,&X1,&Y1); X=X1;Y=Y1;n-=2; }
      int col=R-abs(b)+a; yh=(Y>>(n-1-col))&1; }
    int y0=((w&mW)>>R)&1; if(cst==99) yh=y0^((w*2654435761ULL>>13)%97==0); /* negative control: sparse random g */
    if(yh^y0) G[w>>6]|=1ULL<<(w&63); } }
static void sound(){ G=calloc(nwords(),8); u64 NE=1ULL<<(4*D); u64 tot=0;
  for(int q=0;q<(2*R+1)*(2*R+1)+5;q++){ int a=0,b=0,cst=-1;
    if(q<(2*R+1)*(2*R+1)){ a=q/(2*R+1)-R; b=q%(2*R+1)-R; if(abs(a)+abs(b)>R) continue; }
    else if(q<(2*R+1)*(2*R+1)+4){ cst=q-(2*R+1)*(2*R+1); int x=cst>>1,y=cst&1; if(!(x==y && (x^fb(y,y,y))==y)) continue; }
    else cst=99;
    gmap(a,b,cst); u64 viol=0,fires=0;
    #pragma omp parallel for schedule(dynamic,64) reduction(+:viol,fires)
    for(u64 w=0;w<NW;w++) for(u64 e=0;e<NE;e++){ u64 X0,Y0,wid[64]; int ctx; ext(w,e,&X0,&Y0); patch(X0,Y0,wid,&ctx);
      u64 zm=0; for(int p=0;p<NP;p++) if(!getb(G,wid[p+1])) zm|=1ULL<<p;
      if(crit(ctx,zm,wid)){ fires++; if(getb(G,wid[0])) viol++; } }
    if(cst==99) printf("CONTROL (random sparse g, not commuting; violations expected): "); else if(cst>=0) printf("const %d: ",cst); else printf("sigma^%d F^%d: ",a,b);
    printf("g!=0 on %llu windows, patches where learned criterion fires %llu, violations %llu\n",(unsigned long long)popc(G),(unsigned long long)fires,(unsigned long long)viol);
    if(cst!=99) tot+=viol; fflush(stdout); }
  printf("SOUND rule %d R %d D %d total violations %llu\n",RULE,R,D,(unsigned long long)tot); }
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage\n");return 2;}
  if(!strcmp(argv[1],"minsup")){ setup(atoi(argv[2]),1,atoi(argv[3])); minsup(); u64 tot=0; int mx=0;
    for(int c=0;c<NCTX;c++){tot+=NMS[c]; if(NMS[c]>mx)mx=NMS[c];}
    printf("rule %d D %d points %d ctx %d minimal supports total %llu max %d\n",RULE,D,NP,NCTX,(unsigned long long)tot,mx);
    if(argc>4){ FILE*fp=fopen(argv[4],"w"); fprintf(fp,"%d %d\n",NP,NCTX); for(int p=0;p<NP;p++) fprintf(fp,"%d %d\n",PT[p],PI[p]);
      for(int c=0;c<NCTX;c++){ fprintf(fp,"%d",NMS[c]); for(int k=0;k<NMS[c];k++) fprintf(fp," %llx",(unsigned long long)MS[c][k]); fprintf(fp,"\n"); } fclose(fp); }
    return 0; }
  setup(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
  double t0=now(); EQ=!strcmp(argv[1],"dfse"); minsup(); fprintf(stderr,"minsup %.1fs\n",now()-t0);
  if(!strcmp(argv[1],"sound")){ const char*m=argc>5?argv[5]:"full"; MODE=!strcmp(m,"full")?0:(!strcmp(m,"five")?1:(!strcmp(m,"dc")?3:4)); EQ=!strcmp(m,"dfse"); printf("sound mode %s\n",m); sound(); return 0; }
  MODE=!strcmp(argv[1],"full")?0:(!strcmp(argv[1],"five")?1:(!strcmp(argv[1],"dc")?3:(!strcmp(argv[1],"dfs")||!strcmp(argv[1],"dfse")?4:2)));
  if(MODE==2) loadrules(argv[7]);
  if(MODE==3 && D==2){ u64 n=(u64)NCTX<<NP; DCT=calloc(n/8+1,1);
    #pragma omp parallel for schedule(dynamic,16)
    for(int c=0;c<NCTX;c++) for(u64 z=0;z<(1ULL<<NP);z++) if(dc(c,z)){ u64 k=((u64)c<<NP)|z; __atomic_fetch_or(&DCT[k>>3],(unsigned char)(1u<<(k&7)),__ATOMIC_RELAXED);} 
    fprintf(stderr,"DC table built\n"); }
  return closure(argv[5],argc>6&&strcmp(argv[6],"-")?argv[6]:NULL);
}
