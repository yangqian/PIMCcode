#include "ext.h"
#include "zerorange.h"
#include <gsl/gsl_math.h>
//Below, if not noted, the derivative of normalized relative density matrix 
//are for infinite interaction strength.
//derivative of Ref.2 Eq. 23
real zeroenergy1d(real r,real rp,real coss){
  //trapped version
  if(coss<=0) {
    // printf("energy shouldn't happend\n");
    return 0;}
  real t=r *rp*invsinhepsilon;
  if(t<0.1){
    //series expansion
    return coshvsinhepsilon*(-1+0.5*t-pow2(t)/12.+gsl_pow_4(t)/720.);
  }
  real expf=exp(-t);
  return t*expf*coshvsinhepsilon/(1.-expf);
}
#define Pi 3.141592653589793238
#define LD real
//derivative of Ref.2 Eq. 15 finite g.
real zeroenergy1dfreespacefiniteg(real r,real rp,real coss){
  LD u=(r+rp+g*epsilon)/sqrt(4*epsilon);
  if(coss<=0)
  return
    (g*(2*sqrt(epsilon)*(-(epsilon*g) + r + rp) + expl(pow2(u))*sqrt(Pi)*
        (epsilon*(2 + epsilon*pow2(g)) - pow2(r +
                                              rp))*erfcl(u)))/(8.*epsilon*sqrt(epsilon))
    /(
          1-sqrt(Pi*0.25*epsilon)*g*erfcl(u)*expl(pow2(u))
     )
    ;
  return
    (g*(2*sqrt(epsilon)*(-(epsilon*g) + r + rp) + expl(pow2(u))*sqrt(Pi)*
        (epsilon*(2 + epsilon*pow2(g)) - pow2(r - rp) )*erfcl(u)))/
    (8.*expl((r*rp)/epsilon)*epsilon*sqrt(epsilon))
    /
    (1-expl(-r*rp*invepsilon)*sqrt(Pi*0.25*epsilon)*g*erfcl(u)*expl(pow2(u))
    );
}
//derivative of Ref.2 Eq. 17
real zeroenergy1dfreespace(real r,real rp,real coss){
  if(coss<=0) return 0;
  real t=r *rp *invepsilon;
  return t*exp(-t)*invepsilon/(1-exp(-t));
}
//Ref.2 Eq. 37
real zeroenergy3d(real r,real rp,real coss){
  //trapped version
  real t=r *rp *invsinhepsilon*0.5;
  return -(1+(1+coss)*t)/(exp((1+coss)*t)*t+1)*coshvsinhepsilon;
}
//derivative of Ref.2 Eq. 33
real zeroenergy3dfreespace(real r,real rp,real coss){
  //  return -coshepsilon*(2+(1+coss)*r*rp/sinhepsilon)/2/(exp(0.5*(1+coss)*r*rp/sinhepsilon)+sinhepsilon);
  real t=r *rp *invepsilon*0.5;
  return -(1+(1+coss)*t)/(exp((1+coss)*t)*t+1)*invepsilon;
}
//estimator used for Ref.2, Eq. 9, note all the linear estimator needs a extra
//factor of N, which is canceled out in the calc.c program.
real twob2Z(real x[MAXPARTICLENUM][MAXN][MAXDIMEN]){
  int ii,jj,l,m;
  real old,tempd,newv,newratio=1;
  for(int i=0;i<N;i++){
    old=0;
    newv=0;
    for(ii=0;ii<permute1;ii++){
      for(jj=0;jj<permute2;jj++){
        tempd=1;
        //each term is a product of the single particle terms and pair terms
        //single particle term species 1
        for(l=0;l<identity[0];l++){
          tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][i+1],epsilon);
        }
        //single particle term species 2
        for(l=identity[0];l<ParticleNumber;l++){
          tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][i+1],epsilon);
        }
        newv+=tempd*sig1[ii]*sig2[jj];
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
    newratio*=newv/old;
  }
  return (N*newratio);
}
//twice of trap energy estimator
real PotentialEnergy(real x[MAXLEVYPARTICLENUM][MAXN][MAXDIMEN]){
  int i,j,k;
  real eetemp=0;
  for(i=0;i<ParticleNumber;i++){
    for(j=0;j<N;j+=1){
      for(k=0;k<Dimension;k+=1){
        eetemp+=pow2(x[i][j][k]);
      }
    }
  }
  //eetemp*=2;
  return eetemp;
}
inline real singleparticleenergylink1(real x1[],real x2[]){
  real temp=0;
  for (int k=0;k<Dimension;k++){
    temp-=(pow2(x1[k])+pow2(x2[k]));
  }
  return temp;
}
inline real singleparticleenergylink2(real x1[],real x2[]){
  real temp=0;
  for (int k=0;k<Dimension;k++){
    temp+=x1[k]*x2[k];
  }
  return temp;
}
//thermodynamic energy estimator Ref.3 Eq. 2.109 applied to the pair product
//approximation Re. 3. Eq. 2.34.
real ThermodynamicEnergyInt(real x[MAXPARTICLENUM][MAXN][MAXDIMEN]){
  int ii,jj,l,m;
  real etemp1,etemp2,ene;
  real sumtot=0;
  real nomi,denomi;
  real cp1=0.5*pow2(invsinhepsilon);
  real cp2=(coshepsilon)*pow2(invsinhepsilon);
  real cc=coshepsilon*invsinhepsilon/2*Dimension;
  for(int i=0;i<N;i++){
    nomi=0.;
    denomi=0.;
    for(ii=0;ii<permute1;ii++){
      for(jj=0;jj<permute2;jj++){
        tempd=1;
        //each term contains the single particle terms and pair terms
        //single particle term species 1
        etemp1=0;
        etemp2=0;
        for(l=0;l<identity[0];l++){
          tempd*=trapdensitymatrix(x[l][i],x[permutetable[ii][l]][i+1],epsilon);
          //single particle energy term
          etemp1+=singleparticleenergylink1(x[l][i],x[permutetable[ii][l]][i+1]);
          etemp2+=singleparticleenergylink2(x[l][i],x[permutetable[ii][l]][i+1]);
        }
        //single particle term species 2
        for(l=identity[0];l<ParticleNumber;l++){
          tempd*=trapdensitymatrix(x[l][i],x[permutetable[jj][l]][i+1],epsilon);
          //single particle energy term
          etemp1+=singleparticleenergylink1(x[l][i],x[permutetable[jj][l]][i+1]);
          etemp2+=singleparticleenergylink2(x[l][i],x[permutetable[jj][l]][i+1]);
        }
        //pair terms
        ene=0;
        for(l=0;l<identity[0];l++){
          for(m=identity[0];m<ParticleNumber;m++){
            tempd*=zerorangepairlink4(
                  x[l][i],x[m][i],x[permutetable[ii][l]][i+1],x[permutetable[jj][m]][i+1]);
            ene+=zerorangeenergypairlink4(
                  x[l][i],x[m][i],x[permutetable[ii][l]][i+1],x[permutetable[jj][m]][i+1]);
          }
        }
        denomi+=tempd*sig1[ii]*sig2[jj];
        nomi+=tempd*sig1[ii]*sig2[jj]*(etemp1*cp1+etemp2*cp2+ene);
      }
    }
    sumtot+=nomi/denomi;
  }
  return sumtot+N*ParticleNumber*cc;
}
