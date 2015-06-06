#include <vector>
#include <cstdio> //printf
#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <exception>
#include <cmath>       //pow, exp 
#include <cstdlib>     // atoi, abs
#include <cstring>

//#include <boost/random/mersenne_twister.hpp> //randomness
//#include <boost/random/uniform_real_distribution.hpp> //uniform randomness!

#include "globals.hh"
#include "nucleus.hh"
#include "codex_particle_model.hh"
#include "integrate.hh"
#include "can_be_run.hh"


Proximity_potential::Proximity_potential(double rinSet, double coul1Set, double rcSet, double BSet,double vnConstSet, double muSet){
  m_vnConst=vnConstSet;
  m_mu=muSet;
  m_rin=rinSet;
  m_coul1=coul1Set;
  m_coul2=coul1Set/(2*rcSet);
  m_rc=rcSet;
  m_b=BSet;
}

double Proximity_potential::value(double r) const
{
  double fl=(hbar2/(2*m_mu))*(m_l+1)*m_l;
  //printf("fl:%f\n",fl);
  //printf("l:%i\n",l);
  //printf("mu:%f\n",mu);

  double vcen=fl/pow(r,2);
  double phi;
   
  double vc;
  if(r>m_rc){
    vc=m_coul1/r;
  }
  else{
    vc=m_coul2*(3-pow(r/m_rc,2));
  }
  double zeta=(r-m_rin)/m_b;
  if(zeta<1.2511){
    double help=zeta-2.54;
    phi=-0.5*pow(help,2)-0.0852*pow(help,3);
  }
  else{
    phi=-3.437*exp(-zeta/0.75);
  }
  double vn=m_vnConst*phi;
    
 
  //printf("l internal:%i\n",l);
  //printf("vc:%f\n",vc);
  //printf("vcen:%f\n",vcen);
  //printf("vn:%f\n",vcen);
  //printf("r:%f\n",r);
  return vc+vn+vcen;
}

double Proximity_potential::first_derivative(double r) const
{
  double dvc, dvcen, dvn,dphi;
  dvcen=-hbar2*m_l*(m_l+1)/(m_mu*pow(r,3));
  if(r>m_rc){
    dvc=-m_coul1/pow(r,2);
  }
  else{
    dvc=m_coul2*(-2*r*pow(m_rc,-2));
  }
  double zeta  = (r-m_rin)/m_b;
  if(zeta <= 1.2511){
    double help = zeta - 2.54;
    dphi=(-2*0.5*help-3*0.0852*pow(help,2))/m_b;
  }
  else{
    dphi=3.437/(m_b*0.75)*exp(-zeta/0.75);
  }
  dvn=m_vnConst*dphi;
  return dvc+dvn+dvcen;
}

double Proximity_potential::second_derivative(double r) const
{
  double d2vc,d2vcen,d2vn;
  double zeta = (r-m_rin)/m_b;
  double d2phi;
    
  d2vcen=3*hbar2*(m_l+1)*m_l/(3*m_mu*pow(r,4));

  if(r>m_rc){
    d2vc=m_coul1/pow(r,3);
  }
  else{
    d2vc=-m_coul1/pow(m_rc,3);
  }

  if(zeta <= 1.2511){
    double help = zeta - 2.54;
    d2phi= - 0.5112 * help - 1.0;
  }
  else{
    d2phi = - 6.11022 * exp( - zeta / 0.75);
  }
  d2vn= (m_vnConst/(m_b*m_b)) * d2phi;
  return d2vcen + d2vc + d2vn;
}

void Proximity_potential::print_to_file(char* filename)
{
  //printf("l:%i\n",l);
  double dr=0.1;
  double r=0.5;
  FILE *fout;
  fout = fopen(filename,"w");
  while(r<20){
    fprintf(fout,"%f\t%f\n",r,value(r));
    r+=dr;
  }    
  fclose(fout);    
}

void Proximity_potential::print_derivative_to_file(char* filename)
{
  double dr=0.1;
  double r=0.5;
  FILE *fout;
  fout = fopen(filename,"w");
  while(r<20){
    fprintf(fout,"%f\t%f\n",r,first_derivative(r));
    r+=dr;
  }    
  fclose(fout);    
}



Proximity_potential Codex_particle_model::potential_properties(Nucleus mother, Nucleus evaporation) const
{
  //get constants to determine the potential that "evaporation" has to tunnel through to evaporate from its "mother".
  //Using the "Proximity potential" option from MyCodex.
  //"Vaz et. al." the input file says next to radius correction
    
  double gamma=0.9517 * (1.0 - 1.78260 * pow((mother.N()-mother.Z())/mother.A(),2) );
    
  double A2=evaporation.A();
  double A1=mother.A()-A2; //mass of daughter after evaporation
  double m1=mass_model.mass(mother);
  double m2=mass_model.mass(evaporation);
  double Z2=evaporation.Z();
  double Z1=mother.Z()-Z2; //charge of daugther
  double radiusCorrection=0.0; //Vaz et. al.? Value in CODEX_input.dat, in any case... 0.9 fm used in codex
  double b=1; //[fm], 1 for approximatively all nuclei. 

  double r1=1.28 * pow(A1,1.0/3.0) - 0.76 + 0.8/pow(A1,1.0/3.0)  + radiusCorrection;
  double r2=1.28 * pow(A2,1.0/3.0) - 0.76 + 0.8/pow(A2,1.0/3.0)  + radiusCorrection;

  double rc1=r1 * (1.0 - pow(b / r1,2) );
  double rc2=r2 * (1.0 - pow(b / r2,2) );
  double rin=rc1+rc2; //sum or radii of both nuclei
  double cmean=(rc1*rc2)/rin;
  double rc=1.25*pow(A1,1.0/3.0);   

  double coul1 = Z1*Z2*e2;

  double mu=m1*m2/(m1+m2); 
 
  double vnConst = 4*M_PI * gamma * cmean*b;

  //printf("r0 set:%f\n",r0);
  //printf("rin set:%f\n",rin);
  //printf("rc set:%f\n",rc);
  //printf("mu set:%f\n",mu);
  //printf("vnConst set:%f\n",vnConst);
  //printf("Coul1 set:%f\n",coul1);
  //printf("Coul2 set:%f\n",coul2);



  return Proximity_potential(rin, coul1, rc, b, vnConst,mu);
    
}

std::vector<xFx> Codex_particle_model::potential_min_max(Proximity_potential const & pot, double x, int numChanges) const
{
  //starts at x and goes upwards, looking for sign changes in derivative.
  //returns values at each sign change. NOTE: not recommended to let numChanges be more than 2, since the function will search forever...
  std::vector<xFx> extremum(numChanges);
  double dx=1e-2;
  double df=pot.first_derivative(x);
  int sign=(df > 0) ? 1 : ((df < 0) ? -1 : 0);
  int newsign;
  int changes=0;
  while(changes<numChanges){
    x+=dx;
    if(x>30){
      no_minima_exception ex;
      throw ex;   
    }
    df=pot.first_derivative(x);
    newsign=(df > 0) ? 1 : ((df < 0) ? -1 : 0);
    if(newsign!=sign){
      extremum[changes].f=pot.value(x);
      extremum[changes].x=x;
      changes++;
      sign=newsign;
      //printf("sign:%i\n",sign);
    }
  }
  return extremum;
}

double Codex_particle_model::tunneling(double E, Proximity_potential pot, xFx min, xFx max) const
{
  //we first find intersection between E and V(r), (zeroes of V(r)-E)
  //by Newton's method and a clever initial guess
  //       ^
  //_____/   \_____ E
  //    /     \____ V
  //  \/ 


  //Find inner crossing point: 
  double rinner;

  double r,V;
  double rold,df,f;
  double rmax=max.x;
  double rmin=min.x;
  double rthird=(rmax-rmin)/3.0;
  double Vavg=0.5*(max.f+min.f);
  double rtol=1e-2;
  //printf("E:%f\n",E);
  //if E higher is than average between potential min and max we use one set of initial values
  if(E>Vavg){
    r=rmax-rthird;
  }
  else{
    r=rmin+rthird; 
  }
  //int i=0;
  do{
    //printf("tries:%i\n",i);
    //i++;
    rold=r;
    f=pot.value(r)-E;
    df=pot.first_derivative(r);
    r=r-f/df;

    //we expect to find solution between rmin and rmax...
    if(r<rmin || r>rmax){
      printf("Sort of failed!\n");
      //if not, we set r2 to rmax and decrease it until
      //it passes the crossing point
      r=rmax;
      do{
	r=r-1e-2;
	V=pot.value(r);
	//printf("r:%f\n",r);
      }while(V>E);
      break;
    }
  }while(fabs(r-rold)>rtol);
  //we now think we have the inner crossing point
  rinner=r;

  //find outer crossing point. We now use diffrent initial values!
  double router;
  r=rmax;
  do{
    rold=r;
    f=pot.value(r)-E;
    df=pot.first_derivative(r);
    r=r-f/df;
    //we expect outer crossing point to be beyond rmax...
    if(r<rmax){
      r=rmax;
      do{
	r=r+1e-2;
	V=pot.value(r);
      }while(V>E);
      break;
    }
  }while(fabs(r-rold)>rtol);
  router=r;

  //We now want to integrate between rinner to router with a suitable method! 
  double mu=pot.mu();
  double Gamow=(2*sqrt(2*mu)/hbar)*sqrt(integrate(pot,E,rinner,router));
  return 1.0/(1.0+exp(2*Gamow));
}

double Codex_particle_model::transmission(Nucleus initial, Nucleus final, double E, int l, Proximity_potential pot) const
{
  //returns T_l(E), the transmission coefficient for going from initial to final nucleus by evaporating a particle with orbital angular momentum l and kinetic energy E. E=E_i-B_\nu-E_f

  if(E<=0){
    return 0;
  }
  //printf("E>0!!!\n");
 


  int Ai=initial.A();
  int Zi=initial.Z();
  int Af=final.A();
  int Zf=final.Z();
  if((Zi==Zf) && (Ai-Af==1) && (l==0)){
    //neutrons for l=0 have no barrier to tunnel through
    //and we treat them differently.
    //transmission from rectangular step with depth 40.
    return 4*sqrt(E+40)*sqrt(E)/pow(sqrt(E)+sqrt(E+40),2);
  }
  Nucleus evaporation(Zi-Zf, Ai-Zi-Af+Zf);
  double mu=pot.mu();
  double trans_coef;

  pot.set_l(l);

  std::vector<xFx> extremum;
  try{
  extremum=potential_min_max(pot,0.1,2);
  }
  catch(no_minima_exception e){
    //printf("%s within searched range, setting T=0.\n",e.what());
    return 0;
  }
  if(E<extremum[0].f){
    //can't decay to this!
    return 0;
  }

  //printf("min:%f \n",extremum[0].f);



  if(E<extremum[1].f){
    //calculate transmission coef ~ probability to tunnel. 
    trans_coef=tunneling(E,pot,extremum[0],extremum[1]);
  }
  else{
    //no tunnling needed!
    double secondDeriv=pot.second_derivative(extremum[1].x);
    double hbo=(2*M_PI/hbar)*sqrt(mu/fabs(secondDeriv));//hbar*omega
    //printf("hbo:%f\n",hbo);
    trans_coef=1.0/( 1.0 +exp(-(E - extremum[1].f) * hbo));
  }

  return trans_coef;
}


double Codex_particle_model::spin(Nucleus n) const
{
  //we just know these in a few cases.
  //for now, only needed for the evaporated particles.
  int A=n.A();
  int Z=n.Z();
  if(A==1){
    return 0.5;
  }
  if(A==2 && Z==1){
    return 1;
  }
  if(A==3 && Z==1){
    return 0.5;
  }
  if(A==3 && Z==2){
    return 0.5;
  }
  if(A==4 && Z==2){
    return 0;
  }
}



//this is where we sum up all the decays due to the processes this model describes. record decay and cummulative probability distribution in Rsum_decay.
void Codex_particle_model::Rsum(Nucleus n, int jmax,int lmax,std::vector<std::pair<double,Decay> > & Rsum_decay,  int & counter) const
{
  double Ei=n.E();
  double Ji=n.J();
  double rhoM=rho(n,Ei);
  int Zm=n.Z();
  int Nm=n.N();
  int Zd,Nd;
  double BE,SE,s,Jf,Ef,rhoD,trans_coef;
  Nucleus daughter, evap;
  
  for(int particle=0;particle<=5;particle++){
    int Zevap, Nevap;
    //switch case to make so that we sum lowest probability first.
    switch(particle){
    case 0 : Zevap=2; Nevap=1; break; //3He
    case 1 : Zevap=1; Nevap=2; break; //t
    case 2 : Zevap=1; Nevap=1; break; //d
    case 3 : Zevap=2; Nevap=2; break; //4He
    case 4 : Zevap=1; Nevap=0; break; //p
    case 5 : Zevap=0; Nevap=1; break; //n
    }
    Zd=Zm-Zevap; //daughter Z
    Nd=Nm-Nevap; //daughter N
    daughter.set_Z(Zd);
    daughter.set_N(Nd);
    //if we do not have a mass for the daughter, we skip it.
    //this effectively assigns probability 0 to the decay to that daughter
    //which probably is fine, since it should be very short-lived.
    if(!can_be_run(Zd,Nd)){
      continue;
    }
    evap.set_Z(Zevap);
    evap.set_N(Nevap);
    Proximity_potential pot = potential_properties(n, evap);
    s=spin(evap);
    BE=mass_model.excess_mass(daughter)+mass_model.excess_mass(evap)-mass_model.excess_mass(n);
    //printf("Nev,Zev (BE):%i,%i (%e)\n",Nevap,Zevap,BE);

    //sum over all possible final spin J. Note: J=jf/2
    for(int jf=0;jf<jmax;jf++){
      //get lowest energy threshold for decay to take place.
      Jf=(double)jf/2.0;
      daughter.set_J(Jf);
      SE=BE+intrinsic_energy(daughter);
      //sum over all possible evaporation l.
      //note: what is possible depends on intrinsic spin of evap.
      //and how this couple to final spin of daugther
      for(double S=fabs(Jf-s);S<=fabs(Jf+s);S+=1.0){
	//skip iteration if this S does not give integer l.
	//quick and dirty solution!
	if(Ji-S != (int)(Ji-S)){
	  continue;
	}
	Ef=Ei-BE;
	while(Ef>0){
	  Ef=Ef-dE; //we do this first in loop, so we cant decay to i
	  for(int l=abs(Ji-S);l<=abs(Ji+S);l++){
	    //should probably calculate lmax based on Jmax
	    //this 'if' will have to do until then...
	    if(l>lmax){
	      break;
	    }
	    rhoD=rho(daughter,Ef);
	    trans_coef=transmission(n,daughter,Ei-Ef-BE,l,pot); //Should it be SE rather than BE?
	    //printf("trans_coef: %e\n",trans_coef);
	    if(trans_coef==0){
	      break;
	      //subsequent l will be worse, since higher centrifugal barrier
	    }
	    counter++;
	    //store cummulative R in first and decay info in second.
	    std::pair<double,Decay> to_store;
	    to_store.first=Rsum_decay[counter-1].first+(rhoD/rhoM)*trans_coef;
	    //printf("rhoD/rhoM T: %e\n",(rhoD/rhoM)*trans_coef);
	    //printf("Rsum: %e\n",to_store.first);
	    
            to_store.second.E=Ef;
	    to_store.second.Z=Zevap;
	    to_store.second.N=Nevap;
	    to_store.second.j=jf;
	    to_store.second.l=l;
	    Rsum_decay.push_back(to_store);
	    //if((particle==4 || particle==5 )&& counter%10==0){
	    //  printf("l, counter: %i, %i\n",l,counter);
	    //}
	  }
	}
      }
      //printf("jf:%i\n", jf);
    }
    //printf("Rsum Zev=%i, Nev=%i : %e\n",Zevap,Nevap,Rsum_decay[counter].first);
    //printf("counter: %i\n",counter);
    //printf("----------\n");	    
  }
  
}

