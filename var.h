#include"func.h"

real distancetable[MAXPARTICLENUM][MAXPARTICLENUM][MAXN]={{{0}}};
real dilatedistance;

real tempd;
real etemp;
real T;
int permutestatus;
int BLOCKNUM;
int Dimension;
int ParticleNumber;
int LEVYLENGTH;
int LOGHALFLEVYLENGTH;
int LOGHALFPERMUTELENGTH;
int LEVYPARTICLENUM;
int N;
int PREEQUIL;
int COL;
int maxT;
int LOOP;
int sumcurrent[SUMSPECIES]={0};
real store[SUMSPECIES][MAXLEVEL];

real newx[MAXPARTICLENUM][MAXDIMEN];
int position[MAXPARTICLENUM];
int particle[MAXPARTICLENUM];
real levytrial[MAXLEVYPARTICLENUM][MAXLEVYLENGTH+1][MAXDIMEN];
Ran * myran;
long randstatic;
	int *cor=new int[MAXBLOCKNUM];
	real randn[MAXN+2];
	real epsilon;
	real epsilon2;
	real halfepsilon;
  real sinhepsilon;
  real coshepsilon;
  real invsinhepsilon;
  real coshvsinhepsilon;
	real x[MAXPARTICLENUM][MAXN][MAXDIMEN];
  //store dv for fourth order
	real dv[MAXPARTICLENUM][MAXN][MAXDIMEN];
real ndv[MAXPARTICLENUM][MAXN][MAXDIMEN];
real xcopy[MAXPARTICLENUM][MAXN][MAXDIMEN];
real beta;

int initialtimeslice;
int finaltimeslice;
int cyclelength;
int cycleremainlength;
int check;//1 indicates wiggle move wraps around the path;
int cycle[MAXPARTICLENUM+1];
int cycleremain[MAXPARTICLENUM];
real invepsilon;
real invepsilon2;
int identical[MAXPARTICLENUM][MAXPARTICLENUM+1]={{0}};//0 means identical;
int noidentical[MAXPARTICLENUM][MAXPARTICLENUM+1]={{0}};//0 means identical;
int *identity=new int[MAXPARTICLENUM];
  /*Start Estimator Declare*/
int EstimatorNum;
Estimator *estimator[SUMSPECIES];


/*******levy flight paramters****/
int (*samplingmethod) (real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);
int (*centerofmasssamplingmethod) (real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);
real sigmatable[100];//used in levy bisection
real coshtable[100];//used in levy bisection
  int numprocs, rank;//total number of processors and #
  //char processor_name[MPI_MAX_PROCESSOR_NAME];

real PREFACTOR=10;
int energyonlyflag=0;

int dilationflag=0;
real dilationwidth=0;
real dilationwidthangle=0;

real spreadtable[MAXLEVYLENGTH+1];
real squarespreadtable[MAXLEVYLENGTH+1];
real endweight[MAXLEVYLENGTH+1];
real (*zerorangepairlink)(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p);
real (*zeroenergy)(real r,real rp,real coss);

real g=0;
real invg=0;
int totalpermute;
int permute1;
int permute2;
int *permutetable[MAXPARTICLENUM];
int *sig1;
int *sig2;
