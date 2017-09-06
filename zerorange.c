#include "ext.h"
#include "zerorange.h"
//normalized pair density matrix of various versions
//Ref.2 Eq. 37
real zerorangepairlink3d(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p){
  //mass 1, frequency 1.
  real r2,cos;
  r2=rd*rpd;
  cos=dotproduct(r0,r1,r0p,r1p)
    /r2;
  return zerorangeunitary(rd,rpd,cos,sinhepsilon);
}
//Ref.2 Eq. 33
real zerorangepairlink3dfreespace(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p){
  //mass 1
  real r2,cos;
  r2=rd*rpd;
  cos=dotproduct(r0,r1,r0p,r1p)
    /r2;
  return zerorangeunitary(rd,rpd,cos,epsilon);
}
//Ref.2 Eq. 23
real zerorangepairlink1d(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p){
  //mass 1 trapped verion
  real d=dotproduct(r0,r1,r0p,r1p);
  if(d<0) return 0;
  return 1-exp(-d*invsinhepsilon);
}
#define Pi 3.141592653589793238
#define LD real
//Ref.2 Eq. 15 finite g.
real zerorangepairlink1dfreespacefiniteg
(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p){
  real d=dotproduct(r0,r1,r0p,r1p);
  LD u=(rd+rpd+g*epsilon)/sqrt(4*epsilon);
  if(d<0) 
  return 1-sqrt(Pi*0.25*epsilon)*g*erfc(u)*exp(pow2(u));
  return 1-exp(-d*invepsilon)*sqrt(Pi*0.25*epsilon)*g*erfcl(u)*expl(pow2(u));
}
//Ref.2 Eq. 17
real zerorangepairlink1dfreespace(real rd,real rpd,real *r0,real *r1,real *r0p, real *r1p){
  //mass 1
  real d=dotproduct(r0,r1,r0p,r1p);
  if(d<0) return 0;
  return 1-exp(-d*invepsilon);
}
