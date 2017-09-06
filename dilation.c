#include "ext.h"
#include <gsl/gsl_combination.h>
#include <assert.h>
#include <cmath>
#include "dilation.h"
#include "zerorange.h"
inline real squaresum(real *x){
  int i; real temp=0;
  for (i=0;i<Dimension;i++){
    temp+=x[i]*x[i];
  }
  return temp;
}
//Pair distance move proposal distribution: Ref. 1. Chapter 4.3.3.
real Dilation::proposeLine(int i,int j, int k,real
      x[][MAXN][MAXDIMEN],real newj[],real newk[]){
  real rel[MAXDIMEN];
  real norm;
  int l;
  real delta=(myran->doub()-0.5)*width;
  real tempd=1.;
  dilatedistance=delta+singleslicedistance(x[j][i],x[k][i]);
  norm=0;
  for(l=0;l<Dimension;l++){
    rel[l]=(x[j][i][l]-x[k][i][l]);
    norm+=rel[l]*rel[l];
  }
  norm=sqrt(norm);
  for(l=0;l<Dimension;l++){
    newj[l]=0.5*delta*rel[l]/norm+x[j][i][l];
    newk[l]=-0.5*delta*rel[l]/norm+x[k][i][l];
  }
  if(dilatedistance<0){
    dilatedistance=-dilatedistance;
  }
  return tempd;
}
//a variation of the proposal distribution of the pair distance move.
real Dilation::proposeWithAngle(int i,int j, int k,real x[][MAXN][MAXDIMEN],real newj[],real newk[]){
  real rel[MAXDIMEN];
  real centerofmassjk[MAXDIMEN];
  real norm;
  int l;
  real delta=(myran->doub()-0.5)*width;
  real tempd=1.;
  dilatedistance=delta+singleslicedistance(x[j][i],x[k][i]);
  if(dilatedistance<0)
  dilatedistance=-dilatedistance;
  for(l=0;l<Dimension;l++){
    centerofmassjk[l]=0.5*(x[j][i][l]+x[k][i][l]);
    rel[l]=myran->doub()*2.0-1.0;
  }
  while((norm=squaresum(rel))>1){
    for(l=0;l<Dimension;l++){
      rel[l]=myran->doub()*2.0-1.0;
    }
  }
  norm=0.5*dilatedistance/sqrt(norm);
  for(l=0;l<Dimension;l++){
    newj[l]=norm*rel[l]+centerofmassjk[l];
    newk[l]=-norm*rel[l]+centerofmassjk[l];
  }
  return tempd;
}
//pair distance accept/reject according to acceptance distribution for a
//single pair
int Dilation::dilate(int i,int j, int k,real x[][MAXN][MAXDIMEN]){
  real newj[MAXDIMEN];
  real newk[MAXDIMEN];
  (this->*propose)(i,j,k,x,newj,newk);
  real old=0;
  real newv;
  int l,m;
  real tempd;
  int ii,jj;
  real ratio=1;
  //old term first
  //there are permute1*permute2 total terms
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
  real *pl,*pm;
  newv=0;
  for(ii=0;ii<permute1;ii++){
    for(jj=0;jj<permute2;jj++){
      tempd=1;
      //each term needs the single particle terms and pair terms
      //single particle term species 1
      for(l=0;l<identity[0];l++){
        if(l==j){
          tempd*=trapdensitymatrix(newj,x[permutetable[ii][l]][i+1],epsilon);
        }
        else if(l==k){
          tempd*=trapdensitymatrix(newk,x[permutetable[ii][l]][i+1],epsilon);
        }
        else
        tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][i+1],epsilon);
      }
      //single particle term species 2
      for(l=identity[0];l<ParticleNumber;l++){
        if(l==j){
          tempd*=trapdensitymatrix(newj,x[permutetable[jj][l]][i+1],epsilon);
        }
        else if(l==k){
          tempd*=trapdensitymatrix(newk,x[permutetable[jj][l]][i+1],epsilon);
        }
        else
        tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][i+1],epsilon);
      }
      //pair terms
      for(l=0;l<identity[0];l++){
        if(l==j){
          pl=newj;
        }
        else if(l==k){
          pl=newk;
        }
        else pl=x[l][i];
        for(m=identity[0];m<ParticleNumber;m++){
          if(m==j){
            pm=newj;
          }
          else if (m==k){
            pm=newk;
          }
          else pm=x[m][i];
          tempd*=zerorangepairlink4(
                pl,pm,x[permutetable[ii][l]][i+1],x[permutetable[jj][m]][i+1]);
        }
      }
      newv+=tempd*sig1[ii]*sig2[jj];
    }
  }
  ratio*=(newv/old);
  old=0;
  for(ii=0;ii<permute1;ii++){
    for(jj=0;jj<permute2;jj++){
      tempd=1;
      //each term needs the single particle terms and pair terms
      //single particle term species 1
      for(l=0;l<identity[0];l++){
        tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][(i+N-1)%N],epsilon);
      }
      //single particle term species 2
      for(l=identity[0];l<ParticleNumber;l++){
        tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][(i+N-1)%N],epsilon);
      }
      //pair terms
      for(l=0;l<identity[0];l++){
        for(m=identity[0];m<ParticleNumber;m++){
          tempd*=zerorangepairlink4(
                x[l][i],x[m][i],x[permutetable[ii][l]][(i+N-1)%N],x[permutetable[jj][m]][(i+N-1)%N]);
        }
      }
      old+=tempd*sig1[ii]*sig2[jj];
    }
  }
  newv=0;
  for(ii=0;ii<permute1;ii++){
    for(jj=0;jj<permute2;jj++){
      tempd=1;
      //each term needs the single particle terms and pair terms
      //single particle term species 1
      for(l=0;l<identity[0];l++){
        if(l==j){
          tempd*=trapdensitymatrix(newj,x[permutetable[ii][l]][(i+N-1)%N],epsilon);
        }
        else if(l==k){
          tempd*=trapdensitymatrix(newk,x[permutetable[ii][l]][(i+N-1)%N],epsilon);
        }
        else
        tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][(i+N-1)%N],epsilon);
      }
      //single particle term species 2
      for(l=identity[0];l<ParticleNumber;l++){
        if(l==j){
          tempd*=trapdensitymatrix(newj,x[permutetable[jj][l]][(i+N-1)%N],epsilon);
        }
        else if(l==k){
          tempd*=trapdensitymatrix(newk,x[permutetable[jj][l]][(i+N-1)%N],epsilon);
        }
        else
        tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][(i+N-1)%N],epsilon);
      }
      //pair terms
      for(l=0;l<identity[0];l++){
        for(m=identity[0];m<ParticleNumber;m++){
          if(l==j){
            pl=newj;
          }
          else if (l==k){
            pl=newk;
          }
          else pl=x[l][i];
          if(m==j){
            pm=newj;
          }
          else if (m==k){
            pm=newk;
          }
          else pm=x[m][i];
          tempd*=zerorangepairlink4(
                pl,pm,x[permutetable[ii][l]][(i+N-1)%N],x[permutetable[jj][m]][(i+N-1)%N]);
        }
      }
      newv+=tempd*sig1[ii]*sig2[jj];
    }
  }
  //Eq. 62 of Ref.1
  ratio*=(newv/old);
  if(ratio<0) {
    if(-ratio*pow2(singleslicedistance(newj,newk)/singleslicedistance(x[j][i],x[k][i]))<(*myran).doub()) {
      return(0);}
    else{
      permutestatus^=1;
    }
    //printf("error detected,ratio negative in dilationmove\n");
    //exit(0);
  }
  else if(ratio*pow2(singleslicedistance(newj,newk)/singleslicedistance(x[j][i],x[k][i]))<(*myran).doub()) {
    return(0);}
  for(l=0;l<Dimension;l++){
    x[j][i][l]=newj[l];
  }
  for(l=0;l<Dimension;l++){
    x[k][i][l]=newk[l];
  }
  if(i==0){
    for(l=0;l<Dimension;l++){
      x[j][N][l]=newj[l];
    }
    for(l=0;l<Dimension;l++){
      x[k][N][l]=newk[l];
    }
  }
  return 1;

}
//pair distance move loops over all the beads
int Dilation::update(real x[][MAXN][MAXDIMEN]){
  int ss=0;
  //shuffle the pair order
  shuffle();
  int ii,i,j,k;

  total+=N*psize;
  //loops over all pair
  for(int jj=0;jj<psize;jj++){
    ii=randomarray[jj];
    j=combination[ii*depth];
    k=combination[ii*depth+1];
    assert(j<k);
    //loops over all beads
    for(i=0;i<N;i+=2){
      ss+=dilate(i,j,k,x);
    }
    for(i=1;i<N;i+=2){
      ss+=dilate(i,j,k,x);
    }
  }

  success+=ss;
  return 1;
}
//pair distance move initialize
Dilation::Dilation(real w){
  total=0;
  success=0;
  if (ParticleNumber<2) return;
  if(w>0) { 
    if(!rank){
      printf("Choosing dilation move with Line.\n");
    }
    propose=&Dilation::proposeLine;
    width=w;
  }
  else if(w<0) {
    if(!rank){
      printf("Choosing dilation move with Angle.\n");
    }
    propose=&Dilation::proposeWithAngle;
    width=-w;
  }
  depth=2;
  gsl_combination *c;
  size_t i;
  c=gsl_combination_calloc(ParticleNumber,depth);
  psize=identity[0]*identity[1];
  randomarray=new int[psize];
  for(int k=0;k<psize;k++){
    randomarray[k]=k;
  }
  CombiLength=psize*depth;
  combination=new int[psize*depth];
  i=0;
  do{
    if(
          ((real)gsl_combination_data(c)[0]-identity[0]+0.5)*((real)gsl_combination_data(c)[1]-identity[0]+0.5)>0)
    continue;
    for (int jj=0;jj<depth;jj++){
      *(combination+i+jj)=gsl_combination_data(c)[jj];
    }

    i+=depth;
  }
  while (gsl_combination_next (c)==GSL_SUCCESS);
  assert(i==(size_t)psize*depth);
  // i must at the end of the array
  gsl_combination_free(c);  
}
//get a random order of pairs
void Dilation::shuffle(){
  int temp,i,j;
  for(i=psize-1;i>0;i--){
    j=myran->doub()*(i+1);
    temp=randomarray[j];
    randomarray[j]=randomarray[i];
    randomarray[i]=temp;
  }
}

