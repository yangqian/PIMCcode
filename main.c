#include "var.h"
void parametercheck(){
  //wiggle length
  if(LEVYLENGTH>MAXLEVYLENGTH) {
    if(!rank){
      printf("MAXLEVYLENGTH too small, please redefine\n");
    }
    MPI_Finalize();
    exit(0);
  }
  //path length
  if(N+1>MAXN) {
    if(!rank){
      printf("MAXN too small, please redefine\n");
    }
    MPI_Finalize();
    exit(0);
  }
  //particle number check
  if(ParticleNumber>MAXPARTICLENUM) {
    if(!rank){
      printf("MAXN too small, please redefine\n");
    }
    MPI_Finalize();
    exit(0);
  }
}

int main(int argc, char *argv[]){
  //initialize mpi
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(!rank){
    printf("Total processor used:%d\n",numprocs); 
  }
  //initialize parameters
  initialize(argc,argv);

  //check parameters
  parametercheck();
  //calculate
  calc();

  //finalize
  MPI_Finalize();
}
