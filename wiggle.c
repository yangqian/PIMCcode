#include "ext.h"
#include "zerorange.h"
#include <gsl/gsl_permutation.h>
int factorial(int x){
  int r=1;
  while (x>1){
    r*=x;
    x--;
  }
  return r;
}
//Return the signature of the permutation
//Enable DEBUG flag shows all the permutation within one species
int signp(int n, int P[]){
  if (n==1) return 1;
  int p = 0;
  int v[n];
  int j=n;
  while(j--) {
    v[j] = 0;
  }
  j=n;
  while(j--) {
    if(v[j]) ++p;
    else {
      int x = j;
      do 
      {x = P[x]; v[x] = 1;} 
      while (x!=j);
    }
  }
  return -2*(p&1)+1;
}

int permutationinitialize()
{
  int i,j;
  permute1=factorial(identity[0]);
  permute2=factorial(identity[1]);
  sig1=new int[permute1];
  sig2=new int[permute2];
  sig1[0]=1;
  sig2[0]=1;
  totalpermute=permute1*permute2;
  int maxpermute=permute1>permute2?permute1:permute2;
  for(i=0;i<maxpermute;i++){
    permutetable[i]=new int[ParticleNumber];
    for(j=0;j<ParticleNumber;j++){
      permutetable[i][j]=j;
    }
  }
  if(identity[0]){
    gsl_permutation * p = gsl_permutation_alloc (identity[0]);
    gsl_permutation_init (p);
    int *q=new int[permute1];
    j=0;
    do 
    {
      for(i=0;i<identity[0];i++){
        permutetable[j][i]=gsl_permutation_get(p,i);
        q[i]=gsl_permutation_get(p,i);
      }
      sig1[j]=signp(identity[0],q);
      j++;
    }
    while (gsl_permutation_next(p) ==
          GSL_SUCCESS);
    gsl_permutation_free (p);
    delete q;
  }
  if(identity[1]){
    gsl_permutation * p = gsl_permutation_alloc (identity[1]);
    gsl_permutation_init (p);
    int *q=new int[permute1];
    j=0;
    do 
    {
      for(i=0;i<identity[1];i++){
        permutetable[j][identity[0]+i]=identity[0]+gsl_permutation_get(p,i);
        q[i]=gsl_permutation_get(p,i);
      }
      sig2[j]=signp(identity[1],q);
      j++;
    }
    while (gsl_permutation_next(p) ==
          GSL_SUCCESS);
    gsl_permutation_free (p);

    delete q;
  }
#ifdef DEBUG
  if(!rank)
  for(i=0;i<permute1;i++){
    printf("%d,",sig1[i]);
    for(j=0;j<ParticleNumber;j++){
      printf("%d ",permutetable[i][j]);
    }
    printf("\n");
  }
  if(!rank)
  for(i=0;i<permute2;i++){
    printf("%d,",sig2[i]);
    for(j=0;j<ParticleNumber;j++){
      printf("%d ",permutetable[i][j]);
    }
    printf("\n");
  }
#endif

  return 1;
}



//wiggle move: Ref. 1 4.3.2 Ref. 3 2.5.4
int levyfermi(real x[][MAXN][MAXDIMEN],int
      permutemovetimeslicelength){
  int checkcycle[MAXPARTICLENUM]={0};
  real ratio=1,tempd;
  int i,j=0;
  int ii,jj,l,m;
  real old,newv;
  for(i=0;i<cyclelength;i++)
  checkcycle[cycle[i]]=1;
  for(i=0;i<ParticleNumber;i++){
    if(!checkcycle[i]) cycleremain[j++]=i;
  }

  real (*newpath[MAXPARTICLENUM])[MAXDIMEN];
  cycleremainlength=ParticleNumber-cyclelength;
  int k;
  for(i=0;i<cyclelength;i++){
    for(k=0;k<Dimension;k++){
      levytrial[cycle[i]][0][k]=x[cycle[i]][initialtimeslice][k];//starting point copy
      levytrial[cycle[i]][permutemovetimeslicelength][k]=x[cycle[i]][finaltimeslice][k];
    }
    //Propose new path.
    //Instead of using density matrix in free space as in Ref. 1. Eq. 56,
    //we use the density matrix in trap as in Ref.2 Eq. 35
    levyflight(levytrial[cycle[i]],permutemovetimeslicelength);
    for(j=0;j<permutemovetimeslicelength;j++){
      ratio/=trapdensitymatrix(levytrial[cycle[i]][j],levytrial[cycle[i]][j+1],epsilon);
      ratio*=trapdensitymatrix(x[cycle[i]][(j+initialtimeslice)%N],x[cycle[i]][(j+1+initialtimeslice)%N],epsilon);
    }
    newpath[cycle[i]]=levytrial[cycle[i]];
  }
  for(i=0;i<cycleremainlength;i++){
    newpath[cycleremain[i]]=&x[cycleremain[i]][initialtimeslice];
  }
  if(check){
    for(j=0;j+initialtimeslice<N;j++){
      i=(j+initialtimeslice);
      old=0;
      for(ii=0;ii<permute1;ii++){
        for(jj=0;jj<permute2;jj++){
          tempd=1;
          //Each term contains the single particle terms and pair terms
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
      i=j;
      newv=0;
      for(ii=0;ii<permute1;ii++){
        for(jj=0;jj<permute2;jj++){
          tempd=1;
          //each term contains the single particle terms and pair terms
          //single particle term species 1
          for(l=0;l<identity[0];l++){
            tempd*=trapdensitymatrix(newpath[l][i],newpath[permutetable[ii][l]][i+1],epsilon);
          }
          //single particle term species 2
          for(l=identity[0];l<ParticleNumber;l++){
            tempd*=trapdensitymatrix(newpath[l][i],newpath[permutetable[jj][l]][i+1],epsilon);
          }
          //pair terms
          for(l=0;l<identity[0];l++){
            for(m=identity[0];m<ParticleNumber;m++){
              tempd*=zerorangepairlink4(
                    newpath[l][i],newpath[m][i],newpath[permutetable[ii][l]][i+1],newpath[permutetable[jj][m]][i+1]);
            }
          }
          newv+=tempd*sig1[ii]*sig2[jj];
        }
      }
      ratio*=(newv/old);
    }
    for(ii=0;ii<cycleremainlength;ii++){
      newpath[cycleremain[ii]]-=N;
    }
    for(;j<permutemovetimeslicelength;j++){
      old=0;
      i=(j+initialtimeslice)-N;
      for(ii=0;ii<permute1;ii++){
        for(jj=0;jj<permute2;jj++){
          tempd=1;
          //each term contains the single particle terms and pair terms
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
      i=j;
      newv=0;
      for(ii=0;ii<permute1;ii++){
        for(jj=0;jj<permute2;jj++){
          tempd=1;
          //each term contains the single particle terms and pair terms
          //single particle term species 1
          for(l=0;l<identity[0];l++){
            tempd*=trapdensitymatrix(newpath[l][i],newpath[permutetable[ii][l]][i+1],epsilon);
          }
          //single particle term species 2
          for(l=identity[0];l<ParticleNumber;l++){
            tempd*=trapdensitymatrix(newpath[l][i],newpath[permutetable[jj][l]][i+1],epsilon);
          }
          //pair terms
          for(l=0;l<identity[0];l++){
            for(m=identity[0];m<ParticleNumber;m++){
              tempd*=zerorangepairlink4(
                    newpath[l][i],newpath[m][i],newpath[permutetable[ii][l]][i+1],newpath[permutetable[jj][m]][i+1]);
            }
          }
          newv+=tempd*sig1[ii]*sig2[jj];
        }
      }
      ratio*=(newv/old);
    }
  }
  else
  for(j=0;j<permutemovetimeslicelength;j++){
    old=0;
    i=j+initialtimeslice;
    for(ii=0;ii<permute1;ii++){
      for(jj=0;jj<permute2;jj++){
        tempd=1;
        //each term contains the single particle terms and pair terms
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
    newv=0;
    i=j;
    for(ii=0;ii<permute1;ii++){
      for(jj=0;jj<permute2;jj++){
        tempd=1;
        //each term contains the single particle terms and pair terms
        //single particle term species 1
        for(l=0;l<identity[0];l++){
          tempd*=trapdensitymatrix(newpath[l][i],newpath[permutetable[ii][l]][i+1],epsilon);
        }
        //single particle term species 2
        for(l=identity[0];l<ParticleNumber;l++){
          tempd*=trapdensitymatrix(newpath[l][i],newpath[permutetable[jj][l]][i+1],epsilon);
        }
        //pair terms
        for(l=0;l<identity[0];l++){
          for(m=identity[0];m<ParticleNumber;m++){
            tempd*=zerorangepairlink4(
                  newpath[l][i],newpath[m][i],newpath[permutetable[ii][l]][i+1],newpath[permutetable[jj][m]][i+1]);
          }
        }
        newv+=tempd*sig1[ii]*sig2[jj];
      }
    }
    ratio*=(newv/old);
  }
  //acceptance ratio
  if(ratio<0) {
    if(-ratio<(*myran).doub()) {
      return(0);}
    else{
      //only change sign if it is successful
      permutestatus^=1;
    }
  }
  else if(ratio<(*myran).doub()) {
    return(0);}



  for(i=0;i<cyclelength;i++){
    for(j=1;j<permutemovetimeslicelength;j++){
      for(k=0;k<Dimension;k++){
        x[cycle[i]][(j+initialtimeslice)%N][k]=levytrial[cycle[i]][j][k];
      }
    }
  }
  if(check){
    for(i=0;i<cyclelength;i++){
      for(k=0;k<Dimension;k++){
        x[cycle[i]][N][k]=x[cycle[i]][0][k];
      }
    }
  }
  return 1;
}
//initialize the parameters used for levy bisection in harmonic trap
int levyinitialize(){
  int l;
  real deltaepsilon=epsilon;
  int lmax=(LOGHALFPERMUTELENGTH>LOGHALFLEVYLENGTH?LOGHALFPERMUTELENGTH:LOGHALFLEVYLENGTH);
  real omega=1.0;
  for(l=0;l<=lmax;l++){
    sigmatable[l]=sqrt(0.5*tanh(deltaepsilon * omega));
    coshtable[l]=0.5/cosh(deltaepsilon * omega);
    deltaepsilon*=2;
  }
  real epsilonp;
  real ups1;
  for(l=1;l<LEVYLENGTH;l++){
    epsilonp=epsilon*l;
    ups1=1.0/tanh(epsilon)+1.0/tanh(epsilonp);
    spreadtable[l]=1.0/sqrt(ups1);
    squarespreadtable[l]=1.0/ups1;
    endweight[l]=1.0/sinh(epsilonp);
  }
  return(1);
}
//pick up a particle index and bead index.
int levybisection(real x[][MAXN][MAXDIMEN]){
  int iiii=(int)((*myran).doub()*N);
  cycle[0]=(int)((*myran).doub()*ParticleNumber);
  cyclelength=1;
  initialtimeslice=iiii;
  while(initialtimeslice>=N) initialtimeslice-=N;
  finaltimeslice=initialtimeslice+LEVYLENGTH;
  if(finaltimeslice>=N) {
    check=1;
    finaltimeslice-=N;
  }
  else{
    check=0;
  }
  return levyfermi(x,LEVYLENGTH);
}
//levy flight in harmonic trap 
//levy harmonic path adapted from 
//W. Krauth, Statistical Mechanics: Algorithms and
//Computations, Oxford Master Series in Physics (Oxford University Press,
//Oxford, UK, 2006).
void levyflight(real a[][MAXDIMEN],int n){
  int i,j;
  for (j=0;j<Dimension;j++){
    gaussiansampling(randn,n+1);
    for (i=1;i<n;i++){
      a[i][j]=(a[n][j]*endweight[n-i]+a[i-1][j]*endweight[1])*squarespreadtable[n-i]+randn[i]*spreadtable[n-i];
    }
  }
}
real trapdensitymatrix(real * initial,real * final,real tau){
  //not normalized and lambda=0.5 which means mass is 1 here
  int i;
  real distancesquare=0;
  real dotprod=0;
  for(i=0;i<Dimension;i++){
    dotprod+=initial[i]*final[i];
    distancesquare+=pow2(initial[i]);
    distancesquare+=pow2(final[i]);
  }
  return(exp(-0.5*(distancesquare*cosh(tau)-2*dotprod)/sinh(tau)));
}
