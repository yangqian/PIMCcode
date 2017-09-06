#ifndef __ext__
#define __ext__
#include "func.h"
extern real distancetable[MAXPARTICLENUM][MAXPARTICLENUM][MAXN];
extern real dilatedistance;
extern real tempd;
extern int permutestatus;
extern real T;
extern int BLOCKNUM;
extern int Dimension;
extern int ParticleNumber;
extern int LEVYLENGTH;
extern int LOGHALFLEVYLENGTH;
extern int LOGHALFPERMUTELENGTH;
extern int LEVYPARTICLENUM;
extern int N;
extern int PREEQUIL;
extern int COL;
extern int maxT;
extern int LOOP;
extern int sumcurrent[SUMSPECIES];
extern real store[SUMSPECIES][MAXLEVEL];
extern real newx[MAXPARTICLENUM][MAXDIMEN];
extern int position[MAXPARTICLENUM];
extern int particle[MAXPARTICLENUM];
extern real levytrial[MAXLEVYPARTICLENUM][MAXLEVYLENGTH+1][MAXDIMEN];
extern Ran * myran;
extern long randstatic;
extern int cor[MAXBLOCKNUM];
extern real randn[MAXN+2];
extern real epsilon;
extern real halfepsilon;
extern real epsilon2;
extern real sinhepsilon;
extern real coshepsilon;
extern real invsinhepsilon;
extern real coshvsinhepsilon;
extern real x[MAXPARTICLENUM][MAXN][MAXDIMEN];
//store dv for fourth order
extern real dv[MAXPARTICLENUM][MAXN][MAXDIMEN];
extern real ndv[MAXPARTICLENUM][MAXN][MAXDIMEN];
extern real xcopy[MAXPARTICLENUM][MAXN][MAXDIMEN];
extern real beta;

extern int initialtimeslice;
extern int finaltimeslice;
extern int cyclelength;
extern int cycleremainlength;
extern int check;//1 indicates wiggle move wraps around the path;
extern int cycle[MAXPARTICLENUM+1];
extern int cycleremain[MAXPARTICLENUM];
extern real invepsilon;
extern real invepsilon2;
extern int identical[MAXPARTICLENUM][MAXPARTICLENUM+1];
extern int noidentical[MAXPARTICLENUM][MAXPARTICLENUM+1];
extern int *identity;
/*Start Estimator Declare*/
extern int EstimatorNum;
extern Estimator *estimator[SUMSPECIES];


/*******levy flight paramters****/
extern int (*samplingmethod) (real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);
extern int (*centerofmasssamplingmethod) (real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);
extern real sigmatable[100];
extern real coshtable[100];
extern int numprocs, rank;
extern real etemp;
extern real PREFACTOR;
extern int energyonlyflag;
extern int dilationflag;
extern real dilationwidth;
extern real dilationwidthangle;
extern real spreadtable[MAXLEVYLENGTH+1];
extern real squarespreadtable[MAXLEVYLENGTH+1];
extern real endweight[MAXLEVYLENGTH+1];
inline void SWAP(real &a, real &b) {
  real tempint;tempint=a;a=b;b=tempint;
}
inline void SWAP(int &a, int &b) {
  int tempint;tempint=a;a=b;b=tempint;
}
inline real fpair(real x[],real y[]){ 
  real temp=1;
  for (int l=0;l<Dimension;l++){
    temp*=fabs(x[l]-y[l]);
  }
  return temp;
}
inline real fsingle(real x[]){ 
  real temp=1;
  for (int l=0;l<Dimension;l++){
    temp*=exp(-x[l]*x[l]/2);
  }
  return temp;
}
extern real (*zerorangepairlink)(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p);
extern real (*zeroenergy)(real r,real rp,real coss);

extern real g;
extern real invg;
extern int totalpermute;
extern int permute1;
extern int permute2;
extern int *permutetable[MAXPARTICLENUM];
extern int *sig1;
extern int *sig2;
#endif
