#include "ext.h"
#include <string.h>
#include <wordexp.h>
#include <argp.h>
#include "moves.h"
extern WholePathMove wholepathmove;
extern SingleSliceMove singleslicemove;
const char *argp_program_version =
"PIMCcode 1.0";
const char *argp_program_bug_address =
"<github.com/yangqian/PIMCcode>";
/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  int argv;/* ARG1 and ARG2 */
  int v0flag;              /* v0flag */
  long currenttime;              /* currenttime*/
}; 
struct arguments arguments={0};
/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
   */
static struct argp_option options[] =
{
  {"seed",   772, "Long int", 0,
    "default current time"},
  {"dimension",   'D', "int", 0,
    "range from 1d to 3d,default 3d"},
  {"shift",   's', "Double", 0,
    "whole path shift"},
  {"displacement",   'S', "Double", 0,
    "single slice displacement width"},
  {"onedg",   773, "Double", 0,
    "interaction strength"},
  {"dilation", 'd', "Double", 0, "dilate/shrink move,\"-\" means with angle"},
  {"centerofmassnon", 776, 0, 0, "without centerofmass sampling"},
  {"unitarity", 'u', 0, 0, "2b scattering length diverge"},
  {"levy", 'L', 0, 0, "levy flight for harmonic trap with w=1"},
  {"temperature",   't', "Double", 0,
    "temperature of the system"},
  {"blocklength",   'l', "Int", 0,
    "Length for each block"},
  {"initialize", 'i', "Double", 0, "factors for initialization"},
  {"beta",   'b', "Double", 0,
    "inverse temperature of the system"},
  {"Nbeads",   'n', "EvenInt", 0,
    "Number of beads, better to be even"},
  {"Nwiggle",   'w', "IntPowof2", 0,
    "Wiggle beads, power of 2"},
  {"col",   'c', "Int", 0,
    "correlation length, trials to skip before measure"},
  {"Nparticle1",   '0', "Int", 0,
    "Number of particle of species 1"},
  {"Nparticle2",   '1', "Int", 0,
    "Number of particle of species 1"},
  {"energyonly",   'e', "1/0", 0,
    "true=1,false=0"},
  {"helpinfo", 'h',0,0,"help message"},
  {0}
};
/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
   */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{

  switch (key)
  {
    case 'h':
      {
        argp_usage(state);
      }
      break;
    case 'L':
      samplingmethod=levybisection;
      break;
    case 'u':
      arguments.v0flag=1;
      break;
    case 't':
      sscanf(arg,"" fFormat "",&beta);
      beta=1./beta;
      break;
    case 's':
      real width;
      sscanf(arg,"" fFormat "",&width);
      wholepathmove.setwidth(width);
      break;
    case 'S':
      sscanf(arg,"" fFormat "",&width);
      singleslicemove.setwidth(width);
      break;
    case 'b':
      sscanf(arg,"" fFormat "",&beta);
      break;
    case 'd':
      if(dilationflag==1){
        dilationflag++;
        sscanf(arg,"" fFormat "",&dilationwidthangle);
      }
      else{
        dilationflag++;
        sscanf(arg,"" fFormat "",&dilationwidth);
      }
      break;
    case 'D':
      sscanf(arg,"%d",&Dimension);
      if (Dimension>3||Dimension<1) printf("Dimension=%d choose error\n",Dimension);
      break;
    case 'e':
      sscanf(arg,"%d",&energyonlyflag);
      break;
    case 'i':
      real fac;
      sscanf(arg,"" fFormat "",&fac);
      PREFACTOR/=fac;
      break;
    case 'n':
      sscanf(arg,"%d",&N);
      break;
    case 'l':
      sscanf(arg,"%d",&maxT);
      break;
    case 'w':
      sscanf(arg,"%d",&LEVYLENGTH);
      break;
    case 'c':
      sscanf(arg,"%d",&COL);
      break;
    case '0':
      sscanf(arg,"%d",&identity[0]);
      ParticleNumber+=identity[0];
      break;
    case '1':
      sscanf(arg,"%d",&identity[1]);
      ParticleNumber+=identity[1];
      break;
    case 776:
      centerofmasssamplingmethod=centerofmassnon;
      break;
    case 773:
      sscanf(arg,"" fFormat "",&g);
      invg=1/g;
      break;

    case 772:
      sscanf(arg,"%ld",&arguments.currenttime);
      break;

    case ARGP_KEY_ARG:
      arguments.argv +=1;
      break;
    case ARGP_KEY_END:
      break;
    default:
      if(!rank){
        return ARGP_ERR_UNKNOWN;
      }
      else{
        return 0;
      }
  }
  return 0;
}
/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
   that we accept.
   */
static char args_doc[] = "necessary args: [b|t]lnw01cS";

/*
   DOC.  Field 4 in ARGP.
   Program documentation.
   */
static char doc[] =
"PIMCcode -- Finite temperature path integral Monte Carlo code for determination of virial coefficient.";

/*
   The ARGP structure itself.
   */
static struct argp argp = {options, parse_opt, args_doc, doc};

/*
   The main function.
   Notice how now the only function call needed to process
   all command-line options and arguments nicely
   is argp_parse.
   */
void initialize(int argc, char *argv[]){
  //initialize, give default values
  int i,j;
  int count=0;
  time_t t1;
  LEVYPARTICLENUM=1;
  //fixed by hand, could be changed
  BLOCKNUM=80;
  COL=8;
  maxT=1000;
  wholepathmove.setwidth(0);
  Dimension=3;
  ParticleNumber=0;
  //entra identity[2]=0 is needed to ensure the ending.
  identity[0]=0;identity[1]=0;identity[2]=0;
  N=512;
  LEVYLENGTH=8;
  LOOP=4000001;
  samplingmethod=levybisection;//multilevelbisection;
  centerofmasssamplingmethod=centerofmass;

  t1=time(0);
  arguments.currenttime=t1;
#ifdef DEBUG
  arguments.currenttime=1403902829;
#endif
  //read from the program arguments
  argp_parse (&argp, argc, argv, 0, 0, NULL);
  if(!rank){
    printf("EPOCH TIME:%ld\n",arguments.currenttime);
    printf("HUMAN TIME:%s\n",ctime(&arguments.currenttime));
  }
  //initialize random number generator
  myran=new Ran(arguments.currenttime,rank*100);
  if (arguments.v0flag==1){
    if(Dimension==3){
      zerorangepairlink=zerorangepairlink3d;
      zeroenergy=zeroenergy3d;
    }
    else if (Dimension==1){
      zerorangepairlink=zerorangepairlink1d;
      zeroenergy=zeroenergy1d;
    }
  }
  else
  {
    if(Dimension==3){
      if(!rank)
      printf("3dfiniteg version haven't been implemented,g infinity instead\n");
      zerorangepairlink=zerorangepairlink3d;
      zeroenergy=zeroenergy3d;
    }
    else if (Dimension==1){
      if(g!=0){
        if(!rank)
        printf("1dfiniteg freespace used,g=" fFormat "\n",g);
        zerorangepairlink=zerorangepairlink1dfreespacefiniteg;
        zeroenergy=zeroenergy1dfreespacefiniteg;
      }
      else{
        if(!rank)
        printf("Please specify onedg\n");
        MPI_Finalize();
        exit(0);
      }
    }
  }


  int Totalv0Num=count;
  int r0loffset=0; 
  /* Start Estimator Initializtion */
  estimator[r0loffset]=new Estimator;
  estimatorinitialize(estimator[r0loffset],"b2");
  estimator[r0loffset++]->evaluate=twob2Z;
  estimator[r0loffset]=new Estimator;
  estimatorinitialize(estimator[r0loffset],"VirialEnergy");
  estimator[r0loffset++]->evaluate=PotentialEnergy;
  estimator[r0loffset]=new Estimator;
  estimatorinitialize(estimator[r0loffset],"ThermoEnergy");
  estimator[r0loffset++]->evaluate=ThermodynamicEnergyInt;
  EstimatorNum=r0loffset+Totalv0Num;
  //LEVY initialize
  i=-1;
  j=LEVYLENGTH;
  while(j!=1){//get Log of L
    j>>=1;
    i++;
  }
  LOGHALFLEVYLENGTH=i;

  //particle permutation and species initialization
  constructidentical(identity); 
  constructnoidentical(identity); 
  if(!rank){
    printf("Total particle number:%d\n",ParticleNumber);
  }

  permutationinitialize();

  //constants involving epsilon
  epsilon=beta/N;
  T=1/beta;
  invepsilon=1.0/epsilon;
  invepsilon2=invepsilon*invepsilon;
  halfepsilon=epsilon*0.5;
  epsilon2=pow2(epsilon);
  sinhepsilon=sinh(epsilon);
  invsinhepsilon=1/sinh(epsilon);
  coshepsilon=cosh(epsilon);
  coshvsinhepsilon=cosh(epsilon)/(sinh(epsilon));

  levyinitialize();
}
int constructidentical(int a[]){
  int i,j,k=0,l,m,n;
  int array[MAXPARTICLENUM]={0};
  int curpos[MAXPARTICLENUM]={0};
  i=0;j=0;

  while(a[i]!=0){
    l=0;
    while(l<a[i]){
      identical[j++][0]=a[i];
      l++;
    }
    i++;
  }
  i=0;
  n=0;
  while(a[i]!=0){
    for(m=0;m<a[i];m++)
    array[m]=n++;
    for(k=n-a[i];k<n;k++){
      for(m=0;m<a[i];m++)
      identical[k][++curpos[k]]=array[m];
    }
    i++;
  }
  return(1);
}
int constructnoidentical(int a[]){
  //a[] should give a sequence of same sepcies such as 3 1;
  //the noidentical[j][0] indicate how many particles are diffeferent from it;
  //rest indicate the particel number;
  int i,j,k=0,l,m,n;
  int array[MAXPARTICLENUM]={0};
  int curpos[MAXPARTICLENUM]={0};
  i=0;j=0;

  while(a[i]!=0){
    l=0;
    while(l<a[i]){
      noidentical[j++][0]=ParticleNumber-a[i];
      l++;
    }
    i++;
  }
  i=0;
  n=0;
  while(a[i]!=0){
    for(m=0;m<a[i];m++)
    array[m]=n++;
    for(k=0;k<ParticleNumber;k++){
      if((k>=n-a[i])&&(k<n)) continue;
      for(m=0;m<a[i];m++)
      noidentical[k][++curpos[k]]=array[m];
    }
    i++;
  }
  return(1);
}
int estimatorinitialize(Estimator *estimator,const char *name){
  strcpy(estimator->name,name);
  estimator->Etotalfermion=0;
  estimator->Etotal2fermion=0;
  return(1);
}
