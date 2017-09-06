#include "ext.h"//use beta N as external variables
//levy harmonic path direct generation adapted from 
//W. Krauth, Statistical Mechanics: Algorithms and
//Computations, Oxford Master Series in Physics (Oxford University Press,
//Oxford, UK, 2006).
void LHP(real x[MAXN][MAXDIMEN]){//Levy Harmonic path

  real gamma1,gamma2,deltar;
  int i,k;
  deltar=beta/N;
  for(i=0;i<Dimension;i++){
    gaussiansampling(randn);

    for(k=1;k<N;k++){
      gamma1=coth(deltar)+coth((N-k)*deltar);
      gamma2=x[k-1][i]/sinh(deltar)+x[N][i]/sinh((N-k)*deltar);

      x[k][i]=gamma2/gamma1+randn[k]/sqrt(gamma1);
    }
  }
}
