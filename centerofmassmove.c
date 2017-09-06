#include "ext.h"
//center of mass move implementation: Ref.3 Chapter 2.5.8 
int centerofmass(real x[][MAXN][MAXDIMEN]){
  int iiii=(int)((*myran).doub()*N);
  int check,initialtimeslice,finaltimeslice;
  {
    initialtimeslice=iiii;
    int i,j,k,ii,pi,l;
    int l2=1;
    real meanx[MAXN][MAXDIMEN]={0};
    real factor=1./sqrt((double)ParticleNumber);
    for(j=0;j<=N;j++){
      for(k=0;k<Dimension;k++){
        meanx[j][k]=0;
      }
    }
    for(i=0;i<ParticleNumber;i++){
      for(j=0;j<=N;j++){
        for(k=0;k<Dimension;k++){
          meanx[j][k]+=x[i][j][k];
        }
      }
    }
    for(j=0;j<=N;j++){
      for(k=0;k<Dimension;k++){
        meanx[j][k]=(meanx[j][k])*factor;
      }
    }
    finaltimeslice=initialtimeslice+LEVYLENGTH;
    if(finaltimeslice>=N) {
      check=1;
      finaltimeslice-=N;
    }
    else{
      check=0;
    }
    for(k=0;k<Dimension;k++){
      levytrial[0][0][k]=meanx[initialtimeslice][k];//starting point copy
      levytrial[0][LEVYLENGTH][k]=meanx[finaltimeslice][k];
    }
    for(l=LOGHALFLEVYLENGTH;l>=0;l--){
      l2<<=1;
      for(k=0;k<Dimension;k++){
        j=0;
        gaussiansampling(randn,l2>>1);
        for(i=1;i<l2;i+=2){
          j++;

          levytrial[0][i<<l][k]=(levytrial[0][(i-1)<<l][k]+levytrial[0][(i+1)<<l][k])*coshtable[l]+sigmatable[l]*randn[j];
        } 
      } 
    }
    if(check){
      for(ii=initialtimeslice+1,i=1;ii<=N;ii++,i++){
        for(k=0;k<Dimension;k++){
          for(pi=0;pi<ParticleNumber;pi++){
            x[pi][ii][k]+=((levytrial[0][i][k]-meanx[ii][k])*factor);
          }
        }
      }
      for(ii=0,i=N-initialtimeslice;ii<finaltimeslice;ii++,i++){
        for(k=0;k<Dimension;k++){
          for(pi=0;pi<ParticleNumber;pi++){
            x[pi][ii][k]+=((levytrial[0][i][k]-meanx[ii][k])*factor);
          }
        }
      }
      for(i=0;i<ParticleNumber;i++){
        for(k=0;k<Dimension;k++)
        x[i][N][k]=x[i][0][k];
      }
    }
    else{

      for(i=1,ii=initialtimeslice+1;i<LEVYLENGTH;i++,ii++){
        for(k=0;k<Dimension;k++){
          //levytrial[0][i][k]-=meanx[ii][k];
          for(pi=0;pi<ParticleNumber;pi++){
            x[pi][ii][k]+=((levytrial[0][i][k]-meanx[ii][k])*factor);
          }
        }
      }
    }
  }
  return(1);
}
//center of mass move disabled
int centerofmassnon(real x[][MAXN][MAXDIMEN]){
  return(1);
}
