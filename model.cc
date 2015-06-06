#include <stdio.h> //printf
#include <string.h> //string
#include <fstream>
#include <iostream>
#include <sstream>

#include "nucleus.hh"
#include "mass_model.hh"
#include "globals.hh"
#include "model.hh"
#include "can_be_run.hh"


double Model::macroscopic_mass(Nucleus n) const
{
  /*
    Returns the mass deviation given by a macroscopic droplet model. Parameters obtained from Moeller, Nix et. al. 1995.
    This is used to deduce the microscopic contribution by subtracting this macroscopic mass from the experimental mass.
    This microscopic energy contribution is then used to get a correction to the macroscopic density of states. 
    NOTE: this includes ground state deformation energy (macro model is spherical here). 
    http://www.sciencedirect.com/science/article/pii/S0092640X85710029/pdf?md5=b016ce9f6c0442bf5045c147ca7f6520&pid=1-s2.0-S0092640X85710029-main.pdf
  */
  double Z=n.Z();
  double N=n.N();
  double A=n.A();
  double I=(N-Z)/A;

  double e2=1.4399764; //elementary charge squared [MeV fm]
  double MH=7.289034; //hydrogen atom mass excess[MeV]
  double Mn=8.07143; //neutron mass excess

  double ael=1.433E-05; //electron-binding constant [MeV]
  double K=240; //nuclear compressibility constant [MeV]
  double rp=0.8; //proton rms radius [fm]
  double r0=1.16; //nuclear-radius constant [fm]
  double a=0.68; //yukawa plus exponential potential range [fm]
  double adem=0.70; //range of yukawa function used generate nuclear chrage distr. [fm]
  double rmac=4.8; //average pairing-gap constant [MeV]
  double h=6.6; //neutron-proton interaction constant [MeV]
  double W=30; //Wigner term [MeV]
  double L=0; //density-assymetry constant [MeV]
  double a3=0; //curvature energy constan [MeV]


  double a1=16.247; //volume energy constant [MeV]
  double a2=22.92; //surface energy constant [MeV]
  double J=32.73; //symmetry-energy constant [MeV]
  double Q=29.21; //effective surface-stiffness constant [MeV]
  double a0=0; //A0 constant
  double A0=0; //since it doesn't matter...
  double ca=0.436; //charge-asymmetry constant [MeV]
  double C=60; //preexponential compressibility-term constant [MeV]
  double gamma=0.831; //exponential compressibility-term range constant
    
  //relevant integrals for spherical nuclei.
  //MAY WANT TO CHANGE TO USE SHAPE PARAMS IN NUCLEUS CLASS
  //BUT ALGORITHM CURRENTLY ASSUMES MICROSCOPICAL GROUND STATE ENERGY 
  //INCLUDES DEFORMANTION CONTRIBUTION, WHICH ARE SUBTRACTED BY HAND 
  //ACCORDING TO USER SPECIFIED GROUND STATE DEFORMATION ENERGY
  double Bs=1; 
  double Bv=1;
  double Bw=1;
  double Bk=1; 
  double Br=1; 
 
  double x0=r0*pow(A,1.0/3.0)/a;
  double y0=r0*pow(A,1.0/3.0)/adem;
  //for spherical:
  double B1= 1-3*pow(x0,-2.0)+(1+x0)*(2+3*pow(x0,-1)+3*pow(x0,-2))*exp(-2*x0);
  double B2=1-(1+2*x0-x0*x0)*exp(-2*x0);
  double B3= 1-5*pow(y0,-2.0)*(1-(15.0/8.0)*pow(y0,-1)+(21.0/8.0)*pow(y0,-3)
			       -(3.0/4.0)*(1+(9.0/2.0)*pow(y0,-1.0)+7*pow(y0,-2.0)+(7.0/2.0)*pow(y0,-3))*exp(-2*y0)); //coulomb diffuseness factor
  double B4=1+5*(-3*pow(y0,-2.0)+(15.0/2.0)*pow(y0,-3.0)-(63.0/4.0)*pow(y0,-5.0)
		 +(3.0/4.0)*(2*pow(y0,-1.0)+12*pow(y0,-2.0)+32*pow(y0,-3.0)+42*pow(y0,-4.0)+21*pow(y0,-5.0)))*exp(-2*y0);
  //derived constants
  double c1=(3.0/5.0)*(e2/r0);
  double c2=(1/336.0)*(1.0/J + 18.0/K)*c1*c1;
  double c4=(5.0/4.0)*pow(3.0/(2*M_PI),2.0/3.0)*c1;
  double c5=(1/(64.0*Q))*c1*c1;
  double f0=(-1.0/8.0)*(145.0/48.0)*e2*rp*rp/pow(r0,3);

  double Wigner;
  //calculate Wigner energy
  if((int)Z==(int)N && (int)Z%2){
      //Z, N odd and equal
      Wigner=W*(fabs(I)+pow(A,-1.0));
  }
  else{ 
    Wigner=W*fabs(I);
  }

  //calculate average pairing
  double pairing;
  if((int)Z%2 && (int)N%2){
    pairing=(rmac*Bs)*pow(N,-1.0/3.0)+(rmac*Bs)*pow(Z,-1.0/3.0)-h/(Bs*pow(A,2.0/3.0));//odd-odd
  }
  else if((int)Z%2 && (int)N%2==0){
    pairing=(rmac*Bs)*pow(Z,-1.0/3.0); //Z odd, N even
  }
  else if((int)Z%2==0 && (int)N%2){
    pairing=(rmac*Bs)*pow(N,-1.0/3.0); //Z even, N odd
  }
  else if((int)Z%2==0 && (int)N%2==0){
    pairing=0; //even-even
  }

    
  double deltaAvg=(I+(3.0/16.0)*(c1/Q)*(Bv*Bs/B1)*Z*pow(A,-2.0/3.0))/(1+(9.0/4.0)*(J/Q)*(Bs*Bs/B1)*pow(A,-1.0/3.0));
  double epsilonAvg=(C*exp(-gamma*pow(A,1.0/3.0))-2*a2*B2*pow(A,-1.0/3.0)+L*pow(deltaAvg,2.0)+c1*B4*pow(Z,2.0)*pow(A,-4.0/3.0))/K;

  double EMoeller=MH*Z+Mn*N 
    + (-a1+J*pow(deltaAvg,2)-0.5*K*pow(epsilonAvg,2))*A
    +(a2*B1+2.25*(pow(J,2)/Q)*pow(Bs*deltaAvg,2)/B1)*pow(A,2.0/3.0)
    +a3*pow(A,1.0/3.0)*Bk
    +a0*A0
    +c1*pow(Z,2)*pow(A,-1.0/3.0)*B3
    -c2*pow(Z,2)*pow(A,1.0/3.0)*Br
    -c4*pow(Z,4.0/3.0)*pow(A,-1.0/3.0)
    -c5*pow(Z,2)*Bw*Bs*pow(B1,-1.0)
    +f0*pow(Z,2)*pow(A,-1.0)
    -ca*(N-Z)
    +Wigner
    +pairing
    -ael*pow(Z,2.39);

  return EMoeller;
}

double Model::gs_deformation_energy(Nucleus n) const
{
  return 0.0; //best model! 
  //may want to use value in nucleus, but it may be uninitiated. Do this for now...
}

double Model::microscopic_mass(Nucleus n) const
{
  //subtract macroscopic (for a sphere) energy, deformation energy of ground states, from experimental value
  //and thus supposedly gets energy contribution from pairing and shell effects
  return mass_model.excess_mass(n)-macroscopic_mass(n)-gs_deformation_energy(n);
}

double Model::gs_beta(Nucleus n) const
{
  double beta2;
  //if nothing specified by the user, we throw a fit since MyCodex had no model to calculate this!
  //TODO: If applicable to given nuclei, could use experimental data. See for example http://www.nndc.bnl.gov/databases/databases.html .
  if((int)n.gs_beta()==-1){
    printf("\nZ=%i, A=%i has no ground state quadrupole beta2 specified, nor does the model know how to calculate it!\n",n.Z(),n.A());
    beta2=-1;
  }
  else
    beta2=n.gs_beta();
  return beta2;
  //returns the shape of the nucleus n. If a shape is specified, uses that.
  //If no shape is specified, calculates one based on THEORY?
}


Primary_inertia Model::primary_inertia(Nucleus n) const
{
  //gets the inertia tensor given a quadrupole deformation \beta_2=\beta_20, i.e. under the assumption of axial and R (rotation by pi) symmetry.
  int A=n.A();
  //according to Hasse and Myers, Macroscopic Nuclear Physics (1988).
  Primary_inertia inertia;
  double alpha2=n.gs_alpha();
  double r0=1.16; //[fm]
  double i0=(2.0/5.0)*pow(A,5.0/3.0)*u*pow(r0,2); //moment of inertia of sphere in hbar^2
  inertia.x=i0*(1-alpha2+(3.0/7.0)*alpha2*alpha2);
  inertia.y=i0*(1+(1.0/2.0)*alpha2+(9.0/7.0)*alpha2*alpha2);
  inertia.z=inertia.y;
  return inertia;
}

double Model::yrast(Nucleus n) const
{
  //according to thesis
  double Erot;
  int J=n.J();
  Primary_inertia i=primary_inertia(n);
  double ipar=i.x;
  double iperp=i.y;
  Erot=hbar2*J*(J+1)/(2*iperp); //rotational perpendicular?
  return Erot;
}

double Model::level_density_parameter(Nucleus n) const
{
  //according to Hasse and Myers, Macroscopic Nuclear Physics (1988).
  int A=n.A();
  double alpha2=n.gs_alpha();
  double Bs=1+(2.0/5.0)*pow(alpha2,2)-(4.0/105.0)*pow(alpha2,3)-(66.0/175.0)*pow(alpha2,4);
  double Bk=1+(2.0/5.0)*pow(alpha2,2)+(16.0/105.0)*pow(alpha2,3)-(82.0/175.0)*pow(alpha2,4);
  double a=0.06845*A*(1+3.114*pow(A,-1.0/3.0)*Bs+5.626*pow(A,-2.0/3.0)*Bk);
  return a;
}

double Model::pairing(Nucleus n) const
{
  //Changes from Codex: changed pairing refrence energy to be zero for even-even, consistent with the macroscopic model. MyCodex seem inconsistent here. 
  //Change 2: made it so that the average can exclude nuclei for which no data exists. There is a risk that it will be bad if not enough nuclei are averaged over, though...
  //AVERAGE SHELL EFFECT INCLUDED IN FRDM. Damp that too? (not done here).
  int N=n.N();
  int Z=n.Z();
  Nucleus plus1=n;
  Nucleus minus1=n;
  double mm=mass_model.excess_mass(n)-macroscopic_mass(n);
  double mp1m=0;
  double mm1m=0;
  int elible_nuclei=0;
  plus1.set_Z(Z+1);
  minus1.set_Z(Z-1);
  if(can_be_run(N,Z+1)){
    mp1m=mass_model.excess_mass(plus1)-macroscopic_mass(plus1)-mm;
    elible_nuclei++;
  }
  if(can_be_run(N,Z-1)){
    mm1m=mass_model.excess_mass(minus1)-macroscopic_mass(minus1)-mm;
    elible_nuclei++;
  }
  
  double dp=elible_nuclei>0 ? (1.0/elible_nuclei)*(mp1m+mm1m) : 0;
  plus1=n;
  minus1=n;
  plus1.set_N(N+1);
  minus1.set_N(N-1);
  mp1m=0;
  mm1m=0;
  elible_nuclei=0;
  if(can_be_run(N+1,Z)){
    mp1m=mass_model.excess_mass(plus1)-macroscopic_mass(plus1)-mm;
    elible_nuclei++;
  }
  if(can_be_run(N-1,Z)){
    mm1m=mass_model.excess_mass(minus1)-macroscopic_mass(minus1)-mm;
    elible_nuclei++;
  }

  double dn=elible_nuclei>0 ? (1.0/elible_nuclei)*(mp1m+mm1m) : 0;
  double pairing;
    
  if((int)Z%2 && (int)N%2){
    pairing=-dn-dp;
  }
  else if((int)Z%2 && (int)N%2==0){
    pairing=-dp;
  }
  else if((int)Z%2==0 && (int)N%2){
    pairing=-dn;
  }
  else if((int)Z%2==0 && (int)N%2==0){
    pairing=0;
  }
  return pairing;
}

double Model::critical_energy(Nucleus n) const
{
  //calculates the critical energy for a given J, as it is done in MyCodex.
  double J=(double)n.J();
  double Ec0=10; //In MeV. Value taken from MyCodex.
  double Jc=12; //in units of hbar. Value taken from MyCodex.
  double Ec;
  if(J<Jc){
    Ec=Ec0*sqrt(1.0 - pow(J/Jc,2));
  }
  else{
    Ec=0;
  }
  return Ec;
}

double Model::shell_damping(Nucleus n, double Eeff) const
{
  //get the shell damping for a given effective excitation energy.
  //according to Schmidt and Morawek (and Codex)
  int A=n.A();
  double a=level_density_parameter(n);
    
  double Ed=0.4*pow(A,4.0/3.0)/a;
  return 1-exp(-Eeff/Ed);
}

double Model::pair_damping(Nucleus n,double Eeff) const
{
  //get pair damping for a given effective excitation energy
  //according to MyCodex
  double Ec=critical_energy(n);
  return Eeff<Ec?1-pow(1-Eeff/Ec,2.0):1;
}


double Model::rho(Nucleus n, double E) const
{
  //get density of states according to Grossjean + Feldmeier
  //interacting Fermi-gas
  int A=n.A();
  double J=n.J();
  double Eyrast=yrast(n);
  //printf("Ey:%e\n",Eyrast);
  double Edef=gs_deformation_energy(n);
  double Eeff=E-Eyrast-Edef; //energy above yrast line
  //Eeff=E;
  if(Eeff<0){
    return 0; //no states below yrast line for given spin! 
  }

  double a=level_density_parameter(n);
  double Ed=0.4*pow(A,4.0/3.0)/a; //MeV
  Primary_inertia I=primary_inertia(n);
  double Ieff=pow(I.z*I.y*I.x,1.0/3.0);
  double shellAndPairing=microscopic_mass(n);
  double pairingE=pairing(n); 
  double shell=shellAndPairing-pairingE; //beware of sign!   
  double shellDamping=shell_damping(n,Eeff);
  double pairDamping=pair_damping(n,Eeff);
  double effectiveShell=shell*shellDamping;
  double effectivePair=pairingE*pairDamping;
  double U=Eeff+effectiveShell+effectivePair; //effective nuclear excitation energy.
  if(U<1e-5){
    U=1e-5; //to prevent rho=nan, which happens for U=0. We take this an approximation of that limit.
  }

  //Grossjean & Feldmeier interation starts here. beta=1/T, where T is the nuclear temperature.
  double aOverBeta=sqrt(a*U); //initial value
  double aOverBeta2;
  for(int i=0;i<20;i++){
    aOverBeta2= a * U * (1.0 - exp(-aOverBeta));
    aOverBeta= sqrt(aOverBeta2);
    //printf("aOverBeta:%f\n",aOverBeta);
  }
  double beta=a/(aOverBeta);
  double sigma=sqrt(Ieff/(hbar2*beta));
  double rho=(1/(sqrt(2)*24))*(1/(sigma*pow(a,1.0/4.0)))*(1/pow(U,5.0/4.0))*exp(2*sqrt(a*U))*((2*J+1)/(2*pow(sigma,2)))*exp(-J*(J+1)/(2*pow(sigma,2)));

  return rho;
}

/*
double Model::rho(Nucleus n, double E) const
{
  //get density of states according to Grossjean + Feldmeier
  //interacting Fermi-gas
  int A=n.A();
  double J=n.J();

  double Eyrast=yrast(n);
  double Edef=gs_deformation_energy(n);
  double Eeff=E-Eyrast-Edef; //energy above yrast line
  if(Eeff<0){
    return 0; //no states below yrast line for given spin! 
  }

  double a=level_density_parameter(n);
  double Ed=0.4*pow(A,4.0/3.0)/a; //MeV
  Primary_inertia I=primary_inertia(n);
  double Ieff=pow(I.z*I.y*I.x,1.0/3.0);
  double shellAndPairing=microscopic_mass(n);
  double pairingE=pairing(n); 
  double shell=shellAndPairing-pairingE; //beware of sign!   
  double shellDamping=shell_damping(n,Eeff);
  double pairDamping=pair_damping(n,Eeff);
  double effectiveShell=shell*shellDamping;
  double effectivePair=pairingE*pairDamping;
  double U=Eeff+effectiveShell+effectivePair; //effective nuclear excitation energy.
  if(U<1e-5){
    U=1e-5; //to prevent rho=nan, which happens for U=0. We take this an approximation of that limit.
  }

  //Grossjean & Feldmeier interation starts here. beta=1/T, where T is the nuclear temperature.
  double aOverBeta=sqrt(a*U); //initial value
  double aOverBeta2;
  for(int i=0;i<20;i++){
    aOverBeta2= a * U * (1.0 - exp(-aOverBeta));
    aOverBeta= sqrt(aOverBeta2);
    //printf("aOverBeta:%f\n",aOverBeta);
  }
  double beta=a/(aOverBeta);
  double S=beta*U+aOverBeta;
  double rho=1/(sqrt(2)*24) * (exp(S)/pow(U,3.0/2.0))*(1/sqrt(beta))*(1-exp(-aOverBeta))/sqrt(1-0.5*U*beta*exp(-aOverBeta)) *(2*J+1)/(Ieff/(hbar2*beta))*exp(-pow(J+0.5,2)/(2*Ieff/(hbar2*beta)));

  return rho;
}
*/


/*
double Model::rho(Nucleus n, double E) const
{
  //get density of states through MyCodex method 7, Feldmeier+ASFAC. 
  int A=n.A();
  double J=n.J();

  double Eyrast=yrast(n);
  double Edef=gs_deformation_energy(n);
  double Eeff=E-Eyrast-Edef; //plus pairing in Schmidt & Morawek? Not in code, though!
  if(Eeff<0){
    return 0; //no states below yrast line for given spin! 
  }

  double a=level_density_parameter(n);
  double Ed=0.4*pow(A,4.0/3.0)/a;
  Primary_inertia I=primary_inertia(n);
  double shellAndPairing=microscopic_mass(n);
  double pairingE=pairing(n); 
  double shell=shellAndPairing-pairingE; //beware of sign!   
  double shellDamping=shell_damping(n,Eeff);
  double pairDamping=pair_damping(n,Eeff);
  double effectiveShell=shell*shellDamping;
  double effectivePair=pairingE*pairDamping;
  double U=Eeff+effectiveShell+effectivePair; //a diffrent effective nuclear excitation energy.
  if(U<1e-5){
    U=1e-5; //to prevent rho=nan, which happens for U=0. We take this an approximation of that limit.
  }

  //Grossjean & Feldmeier interation starts here. beta=1/T, where T is the nuclear temperature.
  double aOverBeta=sqrt(a*U); //initial value
  double aOverBeta2;
  for(int i=0;i<20;i++){
    aOverBeta2= a * U * (1.0 - exp(-aOverBeta));
    aOverBeta= sqrt(aOverBeta2);
    //printf("aOverBeta:%f\n",aOverBeta);
  }
  double beta=1/(aOverBeta/a);
  double S=beta*U+aOverBeta;
  double rho=exp(S)*(1-exp(-aOverBeta))/sqrt(1-0.5*U*beta*exp(-aOverBeta))
    *(1/(sqrt(2)*48))*(2*J+1)*(beta/pow(U,3.0/2.0))*(1/(I.y*sqrt(I.x)));
  double betaFrozen=beta; //so that the 'correction' factor will be 1 by default
  if(Eeff<1){
    //freeze Eff at 1 and calculate dampings for this.
    //this will then be used to get a new beta (betaFrozen here,betaZ in CODEX)
    Eeff=1;
    effectiveShell=shell*shell_damping(n,Eeff);
    effectivePair=pairingE*pair_damping(n,Eeff);
    U=Eeff+effectiveShell+effectivePair;
    if(U<1e-5){
      U=1e-5; //to prevent rho=nan
    }
    //do another iteration
    aOverBeta=sqrt(a*U); //initial value
    aOverBeta2;
    for(int i=0;i<10;i++){
      aOverBeta2= a * U * (1.0 - exp(-aOverBeta));
      aOverBeta= sqrt(aOverBeta2);
    }
    betaFrozen=1/(aOverBeta/a);
  }  
  //printf("rho: %e\n",rho*pow(betaFrozen/beta,3.0/2.0));
  return rho*pow(betaFrozen/beta,3.0/2.0);
}
*/

double Model::intrinsic_energy(Nucleus n) const
{
  //not sure what to call this energy...
  //get energy that the nucleus cannot loose without changing other quantum number.
  //i.e. energy associated with spin, etc. Yrast energy in this case!
  //deformation energy? (=0 since the function always returns that...)
  //And codex doesn't seem to use deformation energy here...
  //return yrast(n)+gs_deformation_energy(n);
  return yrast(n);

}
