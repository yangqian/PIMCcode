#include "ext.h"
#include "zerorange.h"
real dotproduct(real a1[], real a2[], real b1[],real b2[]){
  int i;
  real temp=0;
  for(i=0;i<Dimension;i++){
    temp+=(a2[i]-a1[i])*(b2[i]-b1[i]);
  }
  return(temp);
}

real singleslicedistance(real x0[],real x1[]){
  int k;
  real temp=0;
  for(k=0;k<Dimension;k++){
    temp+=pow2((x1[k]-x0[k]));
  }
  return(sqrt(temp));
}
//calculate a distancetable as well as the permutestatus.
int distanceinitialize(real x[][MAXN][MAXDIMEN])
{
  //initialie distance table
  int i,j,k;
  for(i=0;i<ParticleNumber;i++)
  for(j=i+1;j<ParticleNumber;j++)
  for(k=0;k<=N;k++)//remeber to keep the last node
  distancetable[i][j][k]=singleslicedistance(x[i][k],x[j][k]);
  //calculate the permute status, i.e., the sign of the path
  int ii,jj,l,m;
  real old,tempd,ratio=1;
  for(i=0;i<N;i++){
    old=0;
    for(ii=0;ii<permute1;ii++){
      for(jj=0;jj<permute2;jj++){
        tempd=1;
        //each term needs the single particle terms and pair terms
        //single particle term species 1
        for(l=0;l<identity[0];l++){
          tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][i+1],epsilon);
        }
        //single particle term species 2
        for(l=identity[0];l<ParticleNumber;l++){
          tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][i+1],epsilon);
        }
        //pair terms
        for(l=0;l<identity[0];l++){
          for(m=identity[0];m<ParticleNumber;m++){
            tempd*=zerorangepairlink4(
                  x[l][i],x[m][i],x[permutetable[ii][l]][i+1],x[permutetable[jj][m]][i+1]);
          }
        }
        old+=tempd*sig1[ii]*sig2[jj];
      }
    }
    ratio*=old;
  }
  permutestatus=ratio>=0?0:1;
  return(1);
}
