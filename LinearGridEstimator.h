#include "ext.h"
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics.h>
#include <assert.h>
#include <cmath>
#include "accumulator.h"
inline real gsl_stats_mean_long_double(real data[],int skip,size_t s){
  real a=0;
 for(size_t i=0;i<s;i++){
   a+=data[i];
 }
 return a/s;
}
inline real gsl_stats_sd_long_double(real data[],int skip,size_t s){
  real average=gsl_stats_mean_long_double(data,skip,s);
  real a=0;
 for(size_t i=0;i<s;i++){
   a+=data[i]*data[i];
 }
  return sqrt(a/s-average*average);
}
class DoubGridEstimator: public AccumulateEstimator {
  public:
    DoubGridEstimator
(const int d,const real cosmax=8.,const real cosmin=0,const int den=512);
    real (DoubGridEstimator::*ChoiceofAccuFunc)(int *,int);
    void Accumulator(){(this->*ChoiceofAccumulator)();}
    void finish();
  protected:
    void (DoubGridEstimator::*ChoiceofAccumulator)();
    void FiniteTempAccumulator();
    void ZeroTempAccumulator();
    int depth;
    int density;
    real GridMin;
    real GridMax;
    real slice;
    real *xx;
    real *xp;
    real *xpt;
   // real *x2;
    real *xf;
    real *xpf;
    real *xpft;
    //real *xf2;
    real *position;
    int *combination;
    int psize;
    int CombiLength;
    real HyperRadius(int *l,int NN){
      int i,j;
      real temp=0;
      for(i=0;i<depth-1;i++){
        for(j=i+1;j<depth;j++){
          temp+=pow2(distancetable[l[i]][l[j]][NN]);
        }
      }
      return sqrt(temp/depth);
    }
    real PairDistribution(int *l,int NN){
      return distancetable[l[0]][l[1]][NN];
    }
    int particlecount;
    real masscenter[MAXDIMEN];
    real Density(int *l,int NN){
      real temp=0;
      for(int i=0;i<Dimension;i++)
        temp+=pow2(x[*l][NN][i]);
      return sqrt(temp);
    }
    real RelativeDensity(int *l,int NN){
      int i,j;
      if (particlecount==0){
        for(i=0;i<Dimension;i++){
          masscenter[i]=0;
          for(j=0;j<ParticleNumber;j++){
            masscenter[i]+=x[j][NN][i];
          }
          masscenter[i]/=ParticleNumber;
        }
      }
      particlecount++;
      if (particlecount==ParticleNumber) particlecount=0;
      real temp=0;
      for(i=0;i<Dimension;i++)
        temp+=pow2(x[*l][NN][i]-masscenter[i]);
      return sqrt(temp);
    }
};
class GridEstimator: public AccumulateEstimator {
  public:
    GridEstimator
(const int d,const real cosmax=8.,const real cosmin=0,const int den=512);
    real (GridEstimator::*ChoiceofAccuFunc)(int *,int);
    void Accumulator(){(this->*ChoiceofAccumulator)();}
    void finish();
  protected:
    void (GridEstimator::*ChoiceofAccumulator)();
    void FiniteTempAccumulator();
    void ZeroTempAccumulator();
    int depth;
    int density;
    real GridMin;
    real GridMax;
    real slice;
    long *xx;
    real *xp;
    real *xpt;
   // real *x2;
    long *xf;
    real *xpf;
    real *xpft;
    //real *xf2;
    real *position;
    int *combination;
    int psize;
    int CombiLength;
    real HyperRadius(int *l,int NN){
      int i,j;
      real temp=0;
      for(i=0;i<depth-1;i++){
        for(j=i+1;j<depth;j++){
          temp+=pow2(distancetable[l[i]][l[j]][NN]);
        }
      }
     // if(temp<threeR02) {
     //   printf("accumetelate error" fFormat "\n",sqrt(temp/depth));
     //   exit(0);
     //   }
      return sqrt(temp/depth);
    }
    real PairDistribution(int *l,int NN){
      return distancetable[l[0]][l[1]][NN];
    }
    int particlecount;
    real masscenter[MAXDIMEN];
    real Density(int *l,int NN){
      real temp=0;
      for(int i=0;i<Dimension;i++)
        temp+=pow2(x[*l][NN][i]);
      return sqrt(temp);
    }
    real RelativeDensity(int *l,int NN){
      int i,j;
      if (particlecount==0){
        for(i=0;i<Dimension;i++){
          masscenter[i]=0;
          for(j=0;j<ParticleNumber;j++){
            masscenter[i]+=x[j][NN][i];
          }
          masscenter[i]/=ParticleNumber;
        }
      }
      particlecount++;
      if (particlecount==ParticleNumber) particlecount=0;
      real temp=0;
      for(i=0;i<Dimension;i++)
        temp+=pow2(x[*l][NN][i]-masscenter[i]);
      return sqrt(temp);
    }
};
