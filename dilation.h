#include "ext.h"
#include <gsl/gsl_combination.h>
#include <assert.h>
#include <cmath>
class Dilation{
  public:
    Dilation(real width);
    int update(real x[][MAXN][MAXDIMEN]);
    int dilate(int i,int j,int k,real x[][MAXN][MAXDIMEN]);
    int gsdilate(int i,int j,int k,real x[][MAXN][MAXDIMEN]);
    inline real statistics(){ return (real)success/total;};
  unsigned long long total;
  unsigned long long success;
  protected:
    int depth;//2
    int *combination;
    int psize;//C_ParticleNumber_2
    int CombiLength;
    int *randomarray;
    real width;
  void shuffle();
  real proposeLine(int i,int j,int k,real x[][MAXN][MAXDIMEN],real
  newj[],real newk[]);
  real proposeWithAngle(int i,int j,int k,real x[][MAXN][MAXDIMEN],real
  newj[],real newk[]);
  real (Dilation::*propose)(int,int,int,real x[][MAXN][MAXDIMEN],real
  newj[],real newk[]);
};
