#include "ext.h"
//helps to calculate summation of an array of numbers
int multilevelsuminitialize(int i){
  sumcurrent[i]=0;
  return(1);
}
int multilevelsum(real input,int ith){
  // global variable int current[ith],real store[][MAXLEVEL],
  int i,cur,j;
  real temp;
  sumcurrent[ith]++;
  cur=sumcurrent[ith];
  i=1;
  do{
    if(cur%2){
      temp=input;
      for(j=1;j<i;j++)
      temp+=store[ith][j];
      store[ith][i]=temp;
      return(1);
    }
    i++;
  }while((cur>>=1));
  return(1);
}
int multilevelsuminitialize(){
  int i;
  for(i=0;i<SUMSPECIES;i++){
    multilevelsuminitialize(i);
  }
  return(1);
}
real multilevelsumoutput(int ith){
  int i=1,j;
  real temp=0;;
  j=sumcurrent[ith];
  do{
    if(j%2){
      temp+=store[ith][i];
    }
    i++;
  }while((j>>=1));
  return(temp);
}
