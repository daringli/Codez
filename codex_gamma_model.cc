#include "globals.hh"
#include "nucleus.hh"
#include "codex_gamma_model.hh"
#include <cstdio> //printf
#include <cstdlib>     /* abs */

double Codex_gamma_model::transmission(Nucleus initial, double E, int l) const
{
  //really the strength function times E^(2l+1), see Schmidt and Morawek 1991
  //compare eq 4.1 and 4.2.
  //Works for spherical nuclei giant resonances

  double transmission;
  double Z=initial.Z();
  double N=initial.N();
  double A=Z+N;
  //some empirical values the user may want to be able to change...
  double GammaGiantDipole=5.0;
  double GammaGiantQpole1=5.0;
  double GammaGiantQpole2=5.0;
  double energyGiantDipole=80.0/pow(A,1.0/3.0);
  double energyGiantQpole1=63.0/pow(A,1.0/3.0);
  double energyGiantQpole2=130.0/pow(A,1.0/3.0);
  double DipoleSumRuleProportion=1.0;
  double QpoleSumRuleProportion=1.0;


  double lorentzian;
  if(l==1){
    //E1 contribution to transmission strength
    lorentzian = 1.0/(pow(pow(energyGiantDipole,2) - pow(E,2),2)+pow(E*GammaGiantDipole,2));
    transmission = pow(E,4)*lorentzian*DipoleSumRuleProportion*GammaGiantDipole*(Z*N/A)*(36.2/3.0)/(3.84e6);
  }
  if(l==2){
    //E2 contribution
    double R=1.18*pow(A,1.0/3.0);
    double c1=R*R*0.6*pow(1+2.5/(R*R),2);
    lorentzian =0.5*GammaGiantQpole1/(pow(pow(energyGiantQpole1,2) - pow(E,2),2)+pow(E*GammaGiantQpole1,2)) + 0.5*GammaGiantQpole2/(pow(pow(energyGiantQpole2,2) - pow(E,2),2)+pow(E*GammaGiantQpole2,2));
    transmission = pow(E,6)*lorentzian*QpoleSumRuleProportion*Z*c1*(((e2/pow(hbar,3))/931.5)/5.0)*2.0/M_PI;
  }

  /*
  //just plotting and returning from here on!
  double dE=0.1;
  double Ep=0;
  FILE *fE1, *fE2;
  char filenameE1[]="E1.dat";
  char filenameE2[]="E2.dat";

  fE1 = fopen(filenameE1,"w");
  fE2=fopen(filenameE2,"w");
  while(Ep<50){
    double E1=1.0/(pow(pow(energyGiantDipole,2) - pow(Ep,2),2)+pow(Ep*GammaGiantDipole,2));
    double E2=0.5/(pow(pow(energyGiantQpole1,2) - pow(Ep,2),2)+pow(Ep*GammaGiantQpole1,2)) + 0.5/(pow(pow(energyGiantQpole2,2) - pow(Ep,2),2)+pow(Ep*GammaGiantQpole2,2));
    fprintf(fE1,"%f\t%f\n",Ep,E1);
    fprintf(fE2,"%f\t%f\n",Ep,E2);

    Ep+=dE;
  }    
  fclose(fE1);
  fclose(fE2);
  */
  return transmission;
}

void Codex_gamma_model::Rsum(Nucleus n, int jmax,int lmax,std::vector<std::pair<double,Decay> > & Rsum_decay, int & counter) const
{
  //decay modes we are actually interested in!
  //printf("Nev,Zev:0,0\n");
  double Ei=n.E();
  double Ji=n.J();
  double rhoM=rho(n,Ei);
  double Jf,SE,Ef,rhoD,trans_coef;
  //sum over all possible final spin J. Note: J=jf/2
  for(int jf=0;jf<jmax;jf++){
    //printf("jf: %i\n",jf);
    //get lowest energy threshold for decay to take place.
    Jf=(double)jf/2.0;
    n.set_J(Jf);
    SE=intrinsic_energy(n);
    //quick and dirty way to skip non-integer l
    if(Ji-Jf!=int(Ji-Jf)){
      continue;
    }
    int l_start=abs(Ji-Jf);
    int l_stop=abs(Ji+Jf);
    for(int l=((1<l_start) ? l_start : 1);l<=((l_stop<2) ? l_stop : 2);l++){
      Ef=Ei;
      while(Ef>SE){
	Ef=Ef-dE; //we do this first, so we cant decay to initial state
	rhoD=rho(n,Ef);
	trans_coef=transmission(n,Ei-Ef-SE,l);
	//printf("trans_coef: %f\n",trans_coef);
	counter++;
	//store cummulative R in first and decay info in second.
	std::pair<double,Decay> to_store;
	to_store.first=Rsum_decay[counter-1].first+(rhoD/rhoM)*trans_coef;

	//printf("rhoD/rhoM T: %e\n",(rhoD/rhoM)*trans_coef);
	//printf("Rsum: %e\n",to_store.first);
	to_store.second.E=Ef;
	to_store.second.Z=0;
	to_store.second.N=0;
	to_store.second.j=jf;
	to_store.second.l=l;
	Rsum_decay.push_back(to_store);
      }          		
    }
  }
  //printf("Rsum gamma: %e\n",Rsum_decay[counter].first);
  //printf("counter: %i\n",counter);	  
  //printf("----------\n");
}
