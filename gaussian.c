#include "ext.h"//use N myran as external variables
//Ref.1 Eq. 2.69 2.70
void gaussiansampling(real* array){ //N better be even
  real x1, x2, w;
  int i;
  for(i=1;i<=N;){
    do {
      x1 = 2.0 * (*myran).doub() - 1.0;
      x2 = 2.0 * (*myran).doub() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    array[i++] = x1 * w;
    array[i++] = x2 * w;
  }
}
void gaussiansampling(real* array,int n){ //N better be even
  real x1, x2, w;
  int i;
  for(i=1;i<=n;){
    do {
      x1 = 2.0 * (*myran).doub() - 1.0;
      x2 = 2.0 * (*myran).doub() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    array[i++] = x1 * w;
    array[i++] = x2 * w;
  }

}
