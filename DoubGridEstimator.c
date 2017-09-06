#include "ext.h"
#include <gsl/gsl_combination.h>
#include <gsl/gsl_statistics.h>
#include <assert.h>
#include <cmath>
#include "LinearGridEstimator.h"
//Estimator that store histograms of real values.
DoubGridEstimator::DoubGridEstimator
(const int d,const real max,const real min,const int den)
{
  particlecount=0;
  ChoiceofAccumulator=&DoubGridEstimator::FiniteTempAccumulator;
  density=den;
  depth=d;
  GridMin=min;
  GridMax=max;
  gsl_combination *c;
  size_t i;
  //initialize combinations
  c=gsl_combination_calloc(ParticleNumber,depth);
  psize=1;
  for(int k=0;k<depth;k++){
    psize*=(ParticleNumber-k);
    psize/=(1+k);
  }
  CombiLength=psize*depth;
  combination=new int[psize*depth];
  i=0;
  do{
    //memcpy(combination+i,gsl_combination_data(c),sizeof(int)*depth);
    for (int jj=0;jj<depth;jj++){
      *(combination+i+jj)=gsl_combination_data(c)[jj];
    }

    i+=depth;
  }
  while (gsl_combination_next (c)==GSL_SUCCESS);
  assert(i==(size_t)psize*depth);
  // i must be at the end of the array
  gsl_combination_free(c);  
  //initialize position a.k.a. x-axis bin center
  position=new real[density];
  slice=(max-min)/density;
  for(int j=0;j<density;j++){
    position[j]=min+slice*j+slice/2;
  }
  xx=new real[density*psize]();
  xf=new real[density*psize]();
  //choose corresponding accumulating functions
  if(depth==1){
    if(ParticleNumber>=2)
    ChoiceofAccuFunc=&DoubGridEstimator::RelativeDensity;
    else
    ChoiceofAccuFunc=&DoubGridEstimator::Density;
  }
  else if(depth==2){
    ChoiceofAccuFunc=&DoubGridEstimator::PairDistribution;
  }
  else{
    ChoiceofAccuFunc=&DoubGridEstimator::HyperRadius;
  }
}
void DoubGridEstimator::FiniteTempAccumulator(){
  int i,j,k;
  //k is the pointer to xx, i is the pointer to combination
  for(i=0,k=0;i<CombiLength;i+=depth,k+=density){
    //only accumulate the even beads
    for(int kk=0;kk<N;kk+=2){
      real x=(this->*ChoiceofAccuFunc)(combination+i,kk);
      j=(int)floor((x-GridMin)/slice);
      if(j>=0&&j<density){
        xx[k+j]+=1/pow2(x);
        if(permutestatus%2){
          xf[k+j]-=1/pow2(x);
        }
        else{
          xf[k+j]+=1/pow2(x);
        }
      }
    }
  }
}
//store infomation to a file
void DoubGridEstimator::finish(){
  //k is the pointer to xx, piterator is the iteration counter
  int piterator,k;
  FILE *fpf=NULL;
  char fn[100]="";
  xpf=new real[psize];
  xpft=new real[psize];
  if(!rank)
  {
    sprintf(fn,"fwave%dd.dat",depth);
    if((fpf=fopen(fn,"w"))==NULL){
      printf("cannot open file %s\n",fn);
      exit(0);
    }
    fprintf(fpf,"#x average std");
    for(piterator=0,k=0;piterator<psize;piterator++,k+=density){
      fprintf(fpf," ");
      for(int i=0;i<depth-1;i++){
        fprintf(fpf,"%d_",combination[piterator*depth+i]);
      }
      fprintf(fpf,"%d",combination[(piterator+1)*depth-1]);
    }
    fprintf(fpf,"\n");
  }
  for(int j=0;j<density;j++){
    for(piterator=0,k=0;piterator<psize;piterator++,k+=density){
      xpf[piterator]=(real(xf[k+j]))*2./BLOCKNUM/etemp/N/(GridMax-GridMin)*density/numprocs;
    }
    MPI_Reduce(xpf,xpft,psize,REALMPI,MPI_SUM,0,MPI_COMM_WORLD);
    if(!rank){
      fprintf(fpf,"" eLFormat "\t" eLFormat "\t" eLFormat "",position[j],gsl_stats_mean_long_double(xpft,1,psize),gsl_stats_sd_long_double(xpft,1,psize));
      for(piterator=0,k=0;piterator<psize;piterator++,k+=density){
        fprintf(fpf,"\t" eLFormat "",xpft[piterator]);
      }
      fprintf(fpf,"\n");
    }
  }

  if(!rank){
    fclose(fpf);
  }
  delete xpf;
  delete xpft;
}

