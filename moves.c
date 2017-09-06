#include "ext.h"
#include "zerorange.h"
#include "moves.h"
void SingleSliceMove::finish(){
  real ratio=(real) accepted/total;
  real ratioall;
  MPI_Reduce(&ratio,&ratioall,1,REALMPI,MPI_SUM,0,MPI_COMM_WORLD);
  ratioall/=numprocs;
  if(!rank){
    printf("Singleslice move acceptance ratio=%f\n",ratioall);
  }
}

//single slice move Ref.3. 2.5.3 
int SingleSliceMove::update(){
  total+=1;
  real (*newpath[MAXPARTICLENUM]);
  real old,newv,tempd,ratio=1;
  real newx[MAXDIMEN];
  real shift[MAXDIMEN];
  int i,ii,jj,l,j,m,k;
  int newn=(int)(N*(*myran).doub());
  for (i=0;i<Dimension;i++){
    shift[i]=width*((*myran).doub()-0.5);
  }
  int newpar=(int)(ParticleNumber*(*myran).doub());
  for (i=0;i<ParticleNumber;i++){
    //if(i==newpar) continue;
    newpath[i]=x[i][newn];
  }
  newpath[newpar]=newx;
  for (j=0;j<Dimension;j++){
    newx[j]=x[newpar][newn][j]+shift[j];
  }
  if(newn){
    for(i=newn-1;i<=newn+1;i+=2){
      old=0;
      for(ii=0;ii<permute1;ii++){
        for(jj=0;jj<permute2;jj++){
          tempd=1;
          //each term needs the single particle terms and pair terms
          //single particle term species 1
          for(l=0;l<identity[0];l++){
            tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][newn],epsilon);
          }
          //single particle term species 2
          for(l=identity[0];l<ParticleNumber;l++){
            tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][newn],epsilon);
          }
          //pair terms
          for(l=0;l<identity[0];l++){
            for(m=identity[0];m<ParticleNumber;m++){
              tempd*=zerorangepairlink4(
                    x[l][i],x[m][i],x[permutetable[ii][l]][newn],x[permutetable[jj][m]][newn]);
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
            tempd*=trapdensitymatrix(x[l][i],newpath[permutetable[ii][l]],epsilon);
          }
          //single particle term species 2
          for(l=identity[0];l<ParticleNumber;l++){
            tempd*=trapdensitymatrix(x[l][i],newpath[permutetable[jj][l]],epsilon);
          }
          //pair terms
          for(l=0;l<identity[0];l++){
            for(m=identity[0];m<ParticleNumber;m++){
              tempd*=zerorangepairlink4(
                    x[l][i],x[m][i],newpath[permutetable[ii][l]],newpath[permutetable[jj][m]]);
            }
          }
          newv+=tempd*sig1[ii]*sig2[jj];
        }
      }
      ratio*=(newv/old);
    }
  }
  else{
    for(int count=0,i=1;count<2;count++,i=N-1){
      old=0;
      for(ii=0;ii<permute1;ii++){
        for(jj=0;jj<permute2;jj++){
          tempd=1;
          //each term needs the single particle terms and pair terms
          //single particle term species 1
          for(l=0;l<identity[0];l++){
            tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][newn],epsilon);
          }
          //single particle term species 2
          for(l=identity[0];l<ParticleNumber;l++){
            tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][newn],epsilon);
          }
          //pair terms
          for(l=0;l<identity[0];l++){
            for(m=identity[0];m<ParticleNumber;m++){
              tempd*=zerorangepairlink4(
                    x[l][i],x[m][i],x[permutetable[ii][l]][newn],x[permutetable[jj][m]][newn]);
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
            tempd*=trapdensitymatrix(x[l][i],newpath[permutetable[ii][l]],epsilon);
          }
          //single particle term species 2
          for(l=identity[0];l<ParticleNumber;l++){
            tempd*=trapdensitymatrix(x[l][i],newpath[permutetable[jj][l]],epsilon);
          }
          //pair terms
          for(l=0;l<identity[0];l++){
            for(m=identity[0];m<ParticleNumber;m++){
              tempd*=zerorangepairlink4(
                    x[l][i],x[m][i],newpath[permutetable[ii][l]],newpath[permutetable[jj][m]]);
            }
          }
          newv+=tempd*sig1[ii]*sig2[jj];
        }
      }
      ratio*=(newv/old);
    }
  }
  //Ref. Eq. 2.80
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
  for(k=0;k<Dimension;k++){
    x[newpar][newn][k]=newx[k];
  }
  if(newn==0){
    for(k=0;k<Dimension;k++){
      x[newpar][N][k]=x[newpar][0][k];
    }
  }
  accepted+=1;
  return 1;
}

void WholePathMove::finish(){
  real ratio=(real) accepted/total;
  real ratioall;
  MPI_Reduce(&ratio,&ratioall,1,REALMPI,MPI_SUM,0,MPI_COMM_WORLD);
  ratioall/=numprocs;
  if(!rank){
    printf("Wholepath move acceptance ratio=%f\n",ratioall);
  }
}

//whole path move Ref.3. 2.5.7
int WholePathMove::update(){
  total+=1;
  real (*newpath[MAXPARTICLENUM])[MAXDIMEN];
  real old,newv,tempd,ratio=1;
  real newx[MAXN][MAXDIMEN];
  real shift[MAXDIMEN];
  int ii,jj,l,j,m,k;
  for (int i=0;i<Dimension;i++){
    shift[i]=width*((*myran).doub()-0.5);
  }
  int newpar=(int)(ParticleNumber*(*myran).doub());
  for (int i=0;i<ParticleNumber;i++){
    if(i==newpar) continue;
    newpath[i]=x[i];
  }
  newpath[newpar]=newx;
  for (int i=0;i<=N;i++){
    for (int j=0;j<Dimension;j++){
      newx[i][j]=x[newpar][i][j]+shift[j];
    }
  }
  for(int i=0;i<N;i++){
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
    newv=0;
    for(ii=0;ii<permute1;ii++){
      for(jj=0;jj<permute2;jj++){
        tempd=1;
        //each term needs the single particle terms and pair terms
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
  for(j=0;j<=N;j++){
    for(k=0;k<Dimension;k++){
      x[newpar][j][k]=newx[j][k];
    }
  }
  accepted+=1;
  return 1;
}
