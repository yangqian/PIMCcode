#include "ext.h"
//includes commonly used zerorange functions
//normalized pair density matrix in trap for the following two functions
//Ref.2 Eq. 37
inline real zerorangeunitary(real r,real rp,real cos,real epsilon){
  real temp=2*epsilon/r/rp;
  return(1+temp*exp(-(1+cos)/temp));///exp(-epsilon*0.25*(r*r+rp*rp));
}
inline real zerorangepairlink4(real *r0,real *r1,real *r0p, real *r1p){
  //mass 1, frequency 1.
  real rd=singleslicedistance(r0,r1);
  real rpd=singleslicedistance(r0p,r1p);
  real r2,cos;
  r2=rd*rpd;
  cos=dotproduct(r0,r1,r0p,r1p)
    /r2;
  return zerorangeunitary(rd,rpd,cos,sinhepsilon);
}
//negative of derivative of reduced pair density matrix 
//with respect to \tau devided by itself
inline real zerorangeenergypairlink4(real *r0,real *r1,real *r0p, real *r1p){
  //mass 1, frequency 1.
  real rd=singleslicedistance(r0,r1);
  real rpd=singleslicedistance(r0p,r1p);
  real r2,cos;
  r2=rd*rpd;
  cos=dotproduct(r0,r1,r0p,r1p)
    /r2;
  return zeroenergy(rd,rpd,cos);
}
