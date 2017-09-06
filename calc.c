#include "ext.h"
#include <unistd.h>
#include "LinearGridEstimator.h"
#include "dilation.h"
#include "moves.h"
WholePathMove wholepathmove;
SingleSliceMove singleslicemove;
int calc(){
  //variables initialize
  int centercol=N/LEVYLENGTH+1;
  time_t t2,t3;
  real td;
  clock_t start,end;
  FILE *fpath=NULL;
  int i,j,k,l;
  long acceptmove=0,totalmove=0;
  real temp[SUMSPECIES];
  long *totalofnum=new long[MAXBLOCKNUM];
  long *totalpermutestatus=new long[MAXBLOCKNUM];
  real multiaveragefermion[SUMSPECIES];
  real multierrorfermion[SUMSPECIES];

  //correlation estimators initialize
  AccumulateEstimator *gridlist[20];
  int gridtotal=0;
  levyinitialize();
  if(!energyonlyflag){
    //GridEstimator: scaled correlation
    //first argument: involved particle number
    //second argument: range 
    //scaled single particle correlation
    gridlist[gridtotal++]=new GridEstimator(1,5);

    //unscaled single particle correlation
    if(Dimension==3) 
    gridlist[gridtotal++]=new DoubGridEstimator(1,5);

    //scaled pair correlation
    if(ParticleNumber>=2){
      if(Dimension==1) 
      gridlist[gridtotal++]=new GridEstimator(2,5);
      else
      gridlist[gridtotal++]=new GridEstimator(2,8);
    }
    if(ParticleNumber>=3){
      //scaled 3-body correlation
      gridlist[gridtotal++]=new GridEstimator(3,8);
      if (ParticleNumber>=4){
        gridlist[gridtotal++]=new GridEstimator(ParticleNumber,8);
      }
    }
  }
  for(int currenttime=0;currenttime<MAXBLOCKNUM;currenttime++){
    totalofnum[currenttime]=0;
    totalpermutestatus[currenttime]=0;
  }
  Dilation dilationmove=Dilation(dilationwidth);
  Dilation dilationmoveangle=Dilation(dilationwidthangle);


  //permutestatus initialization
  for(i=0;i<BLOCKNUM;i++)
  totalpermutestatus[i]=0;

  /**********path initialization************/
  char word[20]="fpath.save";
  int checkfpath=0;
  //count number of saved files
  if(!rank){
    DIR           *d;
    struct dirent *dir;
    d = opendir(".");
    if (d)
    {
      while ((dir = readdir(d)) != NULL)
      {
        if(strstr(dir->d_name,word)){
          checkfpath++;
        }
      }
      closedir(d);
    }
    else{
      printf("List directory failed.\n");
    }
  }
  //broadcast number of saved files
  MPI_Bcast(&checkfpath,1,MPI_INT,0, MPI_COMM_WORLD);
  if(!rank){
    if(checkfpath==0){
      printf("No saved configuration found.\n");
    }
    else if(checkfpath==1){
      printf("There is %d saved configuration.\n",checkfpath);
    }
    else{
      printf("There are %d saved configurations.\n",checkfpath);
    }
  }
  size_t lSize=0;
  //if no saved configurations, try to initialize a path from scratch.
  if(checkfpath==0){
    levyinitialize();

    for(i=0;i<ParticleNumber;i++){
      for(j=0;j<Dimension;j++)
      x[i][0][j]=0,x[i][N][j]=0;
      LHP(x[i]);
      for(j=0;j<Dimension;j++)
      x[i][0][j]=x[i][N/2][j],x[i][N][j]=x[i][N/2][j];
      LHP(x[i]);
    }
    PREEQUIL=maxT*BLOCKNUM/PREFACTOR;
    for(int j=0;j<PREEQUIL;j++){
      for(int ck=0;ck<centercol;ck++){
        centerofmasssamplingmethod(x);
      }

      for(i=0;i<COL;i++){
        samplingmethod(x);
      }
      if(wholepathmove.width){
        for(int j=0;j<wholepathmove.col;j++){
          wholepathmove.update();
        }
      }
      if(singleslicemove.width){
        for(int j=0;j<singleslicemove.col;j++){
          singleslicemove.update();
        }
      }
      //dilation on a straight line
      if(dilationflag)
      dilationmove.update(x);
      //dilation with arbitrary rotation
      if(dilationflag==2)
      dilationmoveangle.update(x);
    }

  }
  else{
    //otherwise, load the path
    //Situation 1: all the path contained in fpath.save....
    if(checkfpath==1){
      int total;
      if(!rank){
        if((fpath=fopen("fpath.save","rb"))==NULL) 
        {
          printf("load fpath error\n");
          exit(0);
        }
        //obtain file size:
        fseek (fpath, 0 , SEEK_END);
        lSize = ftell (fpath);
        rewind (fpath);
        total=lSize/sizeof(real)/(N+1)/MAXDIMEN/ParticleNumber;
        printf("fpath.save contains %d\n",total);
      }
      MPI_Bcast(&total,1,MPI_INT,0, MPI_COMM_WORLD);
      checkfpath=total;
      //if only one copy of path from the saved file, then broadcast to all
      if(total==1){
        if(!rank){
          for ( i=0;i<ParticleNumber;i++){
            if (!fread(&x[i][0][0],(N+1)*MAXDIMEN,sizeof(real),fpath)) break;
          }
          fclose(fpath);
        }
        for ( i=0;i<ParticleNumber;i++){
          MPI_Bcast(&x[i][0][0],(N+1)*MAXDIMEN,REALMPI,0, MPI_COMM_WORLD);
        }
        levyinitialize();
      }
      else{
        //otherwise, send one copy of a path to every processor until it runs
        //out, then repeat the process.
        int dest;
        if(!rank){
          for(int pathRank=1;pathRank<checkfpath&&pathRank<numprocs;pathRank++){
            if(feof(fpath))
            break;
            for ( i=0;i<ParticleNumber;i++){
              if (!fread(&x[i][0][0],(N+1)*MAXDIMEN,sizeof(real),fpath)) break;
            }
            dest=pathRank;
            while(dest<numprocs){
              for (i=0;i<ParticleNumber;i++){
                MPI_Send(&x[i][0][0],(N+1)*MAXDIMEN , REALMPI, dest, i, MPI_COMM_WORLD);
              }
              //printf("rank %d send to %d\n",rank,dest);
              //sleep(3);
              dest+=checkfpath;
            }
          }
          for ( i=0;i<ParticleNumber;i++){
            if (!fread(&x[i][0][0],(N+1)*MAXDIMEN,sizeof(real),fpath)) break;
          }
          fclose(fpath);
          dest=checkfpath;
          while(dest<numprocs){
            for (i=0;i<ParticleNumber;i++){
              MPI_Send(&x[i][0][0],(N+1)*MAXDIMEN , REALMPI, dest, i, MPI_COMM_WORLD);
            }
            //printf("rank %d send to %d\n",rank,dest);
            //sleep(3);
            dest+=checkfpath;
          }
        }
        if(rank){
          for (i=0;i<ParticleNumber;i++){
            MPI_Recv(&x[i][0][0],(N+1)*MAXDIMEN , REALMPI,0,
                  i, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          }
          //printf("rank %d received from %d\n",rank,0);
        }
      }
    }
    //Situation 2: paths contained in individual fpath.save.1 fpath.save.2 ....
    else{
      int dest;
      if(!rank){
        for(int pathRank=1;pathRank<checkfpath&&pathRank<numprocs;pathRank++){
          sprintf(word,"fpath.save.%d",pathRank);
          if((fpath=fopen(word,"rb"))==NULL) 
          {
            printf("load fpath.save.%d error\n",pathRank);
            exit(0);
          }
          for ( i=0;i<ParticleNumber;i++){
            if (!fread(&x[i][0][0],(N+1)*MAXDIMEN,sizeof(real),fpath)) break;
          }
          fclose(fpath);
          dest=pathRank;
          while(dest<numprocs){
            for (i=0;i<ParticleNumber;i++){
              MPI_Send(&x[i][0][0],(N+1)*MAXDIMEN , REALMPI, dest, i, MPI_COMM_WORLD);
            }
            //printf("rank %d send to %d\n",rank,dest);
            dest+=checkfpath;
          }
        }
        if((fpath=fopen("fpath.save","rb"))==NULL) 
        {
          printf("load fpath error\n");
          exit(0);
        }
        //printf("Loading path from fpath.save and send to each node...\n");
        for ( i=0;i<ParticleNumber;i++){
          if (!fread(&x[i][0][0],(N+1)*MAXDIMEN,sizeof(real),fpath)) break;
        }
        fclose(fpath);
        dest=checkfpath;
        while(dest<numprocs){
          for (i=0;i<ParticleNumber;i++){
            MPI_Send(&x[i][0][0],(N+1)*MAXDIMEN , REALMPI, dest, i, MPI_COMM_WORLD);
          }
          //printf("rank %d send to %d\n",rank,dest);
          dest+=checkfpath;
        }
      }
      if(rank){
        for (i=0;i<ParticleNumber;i++){
          MPI_Recv(&x[i][0][0],(N+1)*MAXDIMEN , REALMPI,0,
                i, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        //printf("rank %d received from %d\n",rank,0);
      }
    }
    //printf("communication done on rank %d\n",rank);
    levyinitialize();
    //if all the path are uncorrelated, skip the preequilibration process
  }
  if(checkfpath>=numprocs){
    PREEQUIL=0;
  }
  else
  PREEQUIL=maxT*BLOCKNUM/PREFACTOR;
  //otherwise, start preequilibration.
  start=clock();
  for(int j=0;j<PREEQUIL;j++){
    for(int ck=0;ck<centercol;ck++){
      centerofmasssamplingmethod(x);
    }

    for(i=0;i<COL;i++){
      samplingmethod(x);
    }
    if(wholepathmove.width){
      for(int j=0;j<wholepathmove.col;j++){
        wholepathmove.update();
      }
    }
    if(singleslicemove.width){
      for(int j=0;j<singleslicemove.col;j++){
        singleslicemove.update();
      }
    }
    //dilation on a straight line
    if(dilationflag)
    dilationmove.update(x);
    //dilation with arbitrary rotation
    if(dilationflag==2)
    dilationmoveangle.update(x);
  }
  end=clock();
  t2=time(0);
  if(!rank&&PREEQUIL){
    td=((real)(end-start))*PREFACTOR/CLOCKS_PER_SEC;
    printf("Estimated time for sampling: " fFormat "\n",td);
  }
  //begin calculation
  if(!rank){
    printf("Initial condition: beta=" gFormat "\nN=%d,LEVYL=%d,BlockN=%d,Skipstep=%d,NofaBlock=%d,Dimension=%d\n",beta,N,LEVYLENGTH,BLOCKNUM,COL,maxT,Dimension);
  }
  //write observables to file after one block.
  FILE *colfile[SUMSPECIES];
  char filename[10000]="";
  if(!rank){
    for(int e=0;e<EstimatorNum;e++){
      sprintf(filename,"col%d.dat",e);
      if((colfile[e]=fopen(filename,"w"))==NULL) {
        printf("error in open %s\n",filename);
        exit(0);
      }
    }
  }
  real * realcache= new real[numprocs];
  start=clock();
  for(k=0;k<BLOCKNUM;k++){
    if(k==1&&!rank){
  end=clock();
  if(!rank){
    td=((real)(end-start))*BLOCKNUM/CLOCKS_PER_SEC;
    t3=td+t2;
    printf("Estimated total time: " fFormat " ",td);
    printf("Finish at %s",ctime(&t3));
  }
    }
    multilevelsuminitialize();
    for(i=1;i<=maxT;i++){
      for(int ck=0;ck<centercol;ck++){
        //center of mass move
        centerofmasssamplingmethod(x);
      }
      for(j=0;j<COL;j++){
        //wiggle move
        acceptmove+=samplingmethod(x);
        totalmove+=1;
      }
      if(wholepathmove.width){
        for(j=0;j<wholepathmove.col;j++){
          //whole path move
          wholepathmove.update();
        }
      }
      if(singleslicemove.width){
        for(j=0;j<singleslicemove.col;j++){
          //naive move
          singleslicemove.update();
        }
      }
      //dilation on a straight line
      if(dilationflag)
      dilationmove.update(x);
      //dilation with arbitrary rotation
      if(dilationflag==2)
      dilationmoveangle.update(x);

      //initialize the pair distance table and permutestatus
          distanceinitialize(x);
      /* Accumulate estimator in a block */

      {
        //Accumulate scalor estimators
        for(int e=0;e<EstimatorNum;e++){
          temp[e]=estimator[e]->evaluate(x);
          multilevelsum(temp[e],(e<<1)+2);
        }
        //based on the sign of observable, add or subtract
        if(permutestatus%2){
          totalpermutestatus[k]--;
          for(int e=0;e<EstimatorNum;e++){
            multilevelsum(-temp[e],(e<<1)+3);
          }
        }
        else {
          totalpermutestatus[k]++;
          for(int e=0;e<EstimatorNum;e++){
            multilevelsum(temp[e],(e<<1)+3);
          }
        }
        //Accumulate correlation estimators
        for(int igrid=0;igrid<gridtotal;igrid++){
          gridlist[igrid]->Accumulator();
        }

      }
    }
    //Adds all the estimators in a block up

    for(int e=0;e<EstimatorNum;e++){
      estimator[e]->Efermion[k]=multilevelsumoutput((e<<1)+3);
    }
    //Save scalar estimators to file col0.dat col1.dat...
    for(int e=0;e<EstimatorNum;e++){
      estimator[e]->Efermion[k]=estimator[e]->Efermion[k]/N/maxT;
      MPI_Gather(&(estimator[e]->Efermion[k]),1,REALMPI,realcache,1,REALMPI,0,MPI_COMM_WORLD);
      real etempcache;
      etemp=(real) totalpermutestatus[k]/maxT;
      MPI_Allreduce(&etemp,&etempcache,1,REALMPI,MPI_SUM,MPI_COMM_WORLD);
      if(!rank){
        etemp=etempcache/numprocs;
        fprintf(colfile[e],"" eFormat " ",etemp);
        for(int ee=0;ee<numprocs;ee++){
          fprintf(colfile[e],"" eFormat " ",realcache[ee]);
        }
        fprintf(colfile[e],"\n");
        fflush(colfile[e]);
      }
    }
  }
  //Close files col0.dat col1.dat...
  if(!rank){
    for(int e=0;e<EstimatorNum;e++){
      fclose(colfile[e]);
    }
  }

  etemp=0;
  for(l=0;l<BLOCKNUM;l++)
  etemp+=totalpermutestatus[l];
  etemp/=BLOCKNUM;

  /*Start Estimate Finalize: warning, the error generated in this program
   * is not accurate when number
   * of processor is small. One needs to analyze the raw data from col*.dat*/
  long cccc;
  MPI_Allreduce(&etemp,&temp[0],1,REALMPI,MPI_SUM,MPI_COMM_WORLD);
  etemp=etemp*etemp;
  MPI_Allreduce(&etemp,&temp[1],1,REALMPI,MPI_SUM,MPI_COMM_WORLD);
  MPI_Reduce(&acceptmove,&cccc,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  //temp[2] is sign error for now
  temp[2]=sqrt((temp[1]-temp[0]*temp[0]/numprocs)/numprocs/(numprocs-1));
  if(!rank){
    printf("Permutestatusaveg:" fFormat ",error:" fFormat "\n",temp[0]/numprocs,temp[2]);
  }
  //etemp records the sign
  etemp=temp[0]/numprocs;
  //now records sign percetage error
  temp[2]=temp[2]/etemp;


  for(int e=0;e<EstimatorNum;e++){
    for(l=0;l<BLOCKNUM;l++){
      estimator[e]->Efermion[l]=estimator[e]->Efermion[l]*maxT/etemp;
      estimator[e]->Etotalfermion+=estimator[e]->Efermion[l];
      estimator[e]->Etotal2fermion+=estimator[e]->Efermion[l]*estimator[e]->Efermion[l];

    }
  }

  for(int e=0;e<EstimatorNum;e++){
    estimator[e]->Etotalfermion/=BLOCKNUM;
    estimator[e]->Etotal2fermion=estimator[e]->Etotalfermion*estimator[e]->Etotalfermion;
  }

  for(int e=0;e<EstimatorNum;e++){
    MPI_Reduce(&estimator[e]->Etotalfermion,&multiaveragefermion[e],1,REALMPI,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&estimator[e]->Etotal2fermion,&multierrorfermion[e],1,REALMPI,MPI_SUM,0,MPI_COMM_WORLD);
  }
  if(rank==0){
    real signerror;
    real totalerror;
    for(int e=0;e<EstimatorNum;e++){
      multiaveragefermion[e]/=numprocs;
      multierrorfermion[e]/=numprocs;
      multierrorfermion[e]=sqrt((multierrorfermion[e]-multiaveragefermion[e]*multiaveragefermion[e])/(numprocs-1));
      //sign error
      signerror=multiaveragefermion[e]*temp[2];
      //total error
      totalerror=sqrt(pow2(multierrorfermion[e])+pow2(signerror));
      printf("%s is " gFormat " " gFormat " " gFormat "\n",estimator[e]->name,multiaveragefermion[e],totalerror,multierrorfermion[e]);
    }
  }

  //output statistical information
  if(!rank){
    printf("Wiggle move acceptance ratio=" fFormat "\n",(real)cccc/numprocs/totalmove);
    if(dilationflag)
    printf("Dilate move acceptance ratio=" fFormat ", dilationwidth=" fFormat "\n",dilationmove.statistics(),dilationwidth);
    if(dilationflag==2)
    printf("Dilate move acceptance ratio=" fFormat "a ,dilationwidth=" fFormat "\n",dilationmoveangle.statistics(),dilationwidthangle);
  }
  totalmove=0;acceptmove=0;
  wholepathmove.finish();
  singleslicemove.finish();
  //finalize correlation estimators
  for(int igrid=0;igrid<gridtotal;igrid++){
    gridlist[igrid]->finish();
  }




  MPI_File file;
  MPI_Status status;
  // Opening one shared file to store the path
  MPI_File_open(MPI_COMM_WORLD, "fpath.save",
        MPI_MODE_CREATE|MPI_MODE_WRONLY,
        MPI_INFO_NULL, &file);

  // Setting local offset for each processor
  // Writing the results to the shared file.
  MPI_Offset offset = sizeof(real)*(N+1)*MAXDIMEN*ParticleNumber*rank;
  MPI_File_seek(file, offset, MPI_SEEK_SET);
  for ( i=0;i<ParticleNumber;i++){
    MPI_File_write(file, &x[i][0][0], (N+1)*MAXDIMEN,
          REALMPI, &status);
  }
  MPI_File_close(&file);
  //record the completion time
  if(!rank){
    t3=time(0);
    printf("Finish time: %s",ctime(&t3));
  }

  return(1);

}

