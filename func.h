#ifndef _func_h
#define _func_h
//definition of constants
#define MAXBLOCKNUM 4000
#define MAXMUMOFBETANUM 9 //length of the array store the totalenergy
#define MAXN 32
#define MAXDIMEN 3
#define MAXPARTICLENUM 6//6
#define MAXLEVYLENGTH 4
#define MAXLEVYPARTICLENUM MAXPARTICLENUM
//supports summation of 2^MAXLEVEL times
#define MAXLEVEL 100
//supports at most SUMSPECIES estimators
#define SUMSPECIES 1000
#define VIEWFORMAT "%16.14g\t"
//define ll to choose long double during calculation;nessary for finite g
//calculation.
#ifdef ll
typedef long double real;
#define REALMPI MPI_LONG_DOUBLE
#define gFormat "%19.15Lg"
#define fFormat "%Lf"
#define eFormat "%19.15Le"
#define eLFormat "%19.15Le"
inline long double sqrt(long double x){
  return sqrtl(x);
}
#else
typedef double real;
#define REALMPI MPI_DOUBLE
#define gFormat "%19.15g"
#define fFormat "%lf"
#define eFormat "%19.15e"
#define eLFormat "%19.15e"
#endif
#include <mpi.h>
#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#include<string.h>
//list files in directory
#include <dirent.h> 
#ifndef DEBUG
//turn on production flags for GSL library
#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF
#endif
#include "myran.h"
#define ln(x) log(x)
#define coth(x) (1.0/tanh(x))
//defininion of functions 
extern int calc();
extern void LHP(real x[][MAXDIMEN]);
extern void gaussiansampling(real x[]);
extern void gaussiansampling(real x[],int n);
extern void initialize(int argc, char *argv[]);
extern int permutation(int a[],int);
extern int permutation(int a[]);
extern real freedensitymatrix(real * initial,real * final,real tao);
extern real trapdensitymatrix(real * initial,real * final,real tao);
extern real singleslicedistance(real x0[],real x1[]);
extern int multilevelsum(real input, int ith);
extern real multilevelsumoutput(int ith);
extern int multilevelsuminitialize(int ith);
extern int multilevelsuminitialize();
extern int constructidentical(int a[]);
extern int constructnoidentical(int a[]);
extern int distanceinitialize(real x[][MAXN][MAXDIMEN]);
extern real dotproduct(real a1[], real a2[], real b1[],real b2[]);


typedef struct 
{
  real Efermion[MAXBLOCKNUM];
  real Etotalfermion;
  real Etotal2fermion;
  char name[50];
  real (*evaluate)(real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);//estimator
} Estimator;
extern int estimatorinitialize(Estimator *estimator,const char *name);

extern real PotentialEnergy(real x[MAXLEVYPARTICLENUM][MAXN][MAXDIMEN]);
extern real ThermodynamicEnergyInt(real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);


extern real twob2Z(real x[MAXPARTICLENUM][MAXN][MAXDIMEN]);

extern int levyinitialize();
extern int levybisection(real x[][MAXN][MAXDIMEN]);



extern int centerofmass(real x[][MAXN][MAXDIMEN]);
extern int centerofmassnon(real x[][MAXN][MAXDIMEN]);
extern int dilation(real x[][MAXN][MAXDIMEN]);
inline real pow2(real x){
  return x*x;
}
inline real pow3(real x){
  return x*x*x;
}
void levyflight(real a[][MAXDIMEN],int n);
extern real zeroenergy1d(real r,real rp,real coss);
extern real zeroenergy1dfreespacefiniteg(real r,real rp,real coss);
extern real zeroenergy3d(real r,real rp,real coss);
extern real zerorangepairlink3d(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p);
extern real zerorangepairlink1d(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p);
extern real zerorangepairlink1dfreespacefiniteg(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p);
extern void coordinateupdate(int);

extern int permutationinitialize();
extern int levyfermi(real x[][MAXN][MAXDIMEN],int
      permutemovetimeslicelength);
#endif
