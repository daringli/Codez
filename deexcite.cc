#include <fstream>
#include <vector>


#include <boost/random/uniform_real_distribution.hpp> //uniform randomness!

#include "globals.hh"
#include "deexcite.hh"
#include "decay.hh"
#include "lorentz_boost.hh"

#include "codex_particle_model.hh" //since I can't get B.E. from encapsulated calcs... ;_;

#include <boost/math/special_functions/legendre.hpp>

void deexcite(std::vector<Nucleus> &products, const std::vector<std::pair<double,Decay> > &table)
{
  //given a table of deexcitation probabilities, does one decay of products[0]
  //and adds the resulting evaporated particles to products vector

  long double max=table[table.size()-1].first; //the whole R sum.
  boost::random::uniform_real_distribution<long double> dist(0,max);
  long double random=dist(gen);
  //printf("random (max): %Le (%Le) \n",random, max);
  int i=0;
  while(table[i].first < random){
    i++;
  }
  
  //get Ei out of table => magnitude of momentum.
  double Ei=table[1].second.E; //0 is empty since we need to zero initialize sum
  double Ef=table[i].second.E;
  //printf("Ef: %f\n",Ef);
  
  //This nucleus will update products[0]
  Nucleus daughter;
  daughter.set_N(products[0].N()-table[i].second.N);
  daughter.set_Z(products[0].Z()-table[i].second.Z);
  daughter.set_J(table[i].second.j/2.0);
  daughter.set_E(table[i].second.E);

  Nucleus evap;
  evap.set_N(table[i].second.N);
  evap.set_Z(table[i].second.Z);
  evap.set_E(0.0);


  double p_evap_zm;
  if(evap.N()!=0 && evap.Z()!=0){
    //not photons!
    double mass_evap=mass_model.mass(evap);
  
    //generate momentum of evaporate and daughter (recoil!)
    //evaporation ZM frame, we have
    double gamma_nuc_to_evap = (Ei+ mass_model.mass(products[0])-Ef - mass_model.mass(daughter))/mass_evap;
    double beta_nuc_to_evap = sqrt(1-1/gamma_nuc_to_evap);
    p_evap_zm= gamma_nuc_to_evap*mass_evap*beta_nuc_to_evap;
  }
  else{
    //photons!
    p_evap_zm = Ei+ mass_model.mass(products[0])-Ef - mass_model.mass(daughter);
    //printf("p_gamma: %e\n",p_evap_zm);
  }

  //generate direction from legendre polynomial and l
  int l=table[i].second.l;
  boost::random::uniform_real_distribution<double> dist_phi(0,2*M_PI);
  double phi=dist_phi(gen);
  //integrate P_l^2(x) until we cross this random threshold
  boost::random::uniform_real_distribution<double> urand(0,2.0/(2*l+1));
  double pl2_rng=urand(gen);
  double cost=-1;
  double dcost=0.01;
  double pl2_int=0;
  while(pl2_int<pl2_rng){
    cost=cost+dcost;
    pl2_int=pl2_int+pow(boost::math::legendre_p(l, cost),2);
  }
  double theta=M_PI-acos(cost);

  //assign momentum to evap (we don't return this, since not in right frame yet)
  Vector4 P_evap_zm;
  P_evap_zm.v3.x = p_evap_zm*sin(theta)*cos(phi);
  P_evap_zm.v3.y = p_evap_zm*sin(theta)*sin(phi);
  P_evap_zm.v3.z = p_evap_zm*cos(theta);
  P_evap_zm.x0=sqrt(pow(mass_model.mass(evap),2)+P_evap_zm.v3*P_evap_zm.v3);
  

  Vector4 P_daughter_zm;
  P_daughter_zm.v3=-P_evap_zm.v3;
  P_daughter_zm.x0=sqrt(pow(mass_model.mass(daughter)+Ef,2)+P_evap_zm.v3*P_evap_zm.v3);

    //ZM frame is also the mother frame, which allows us to get relation to lab-frame. beta=p/E. Calculate lab frame momentum.
  Vector3 beta_zm_to_lab=-products[0].P().v3*(1/(products[0].P().x0));
  evap.set_P(lorentz_boost(P_evap_zm,beta_zm_to_lab));
  daughter.set_P(lorentz_boost(P_daughter_zm,beta_zm_to_lab));
  
  //generate random number to spoil E* with to make less discrete
  //boost::random::uniform_real_distribution<double> dist_e(0,2*dE);
  //double random_E=dist_e(gen)-dE; //random number to spoil energy with

  //print info about decay and final state to file, along with momentum of evap.
  //char[] decay_filename="interm/decay_info.tsv";
  //FILE* decay_file=fopen(decay_filename,'w');

  //fprintf(decay_file,"%i\t %i\t %i\t %.1f\t  %.1f\t %.1f\t %.1f\t %.1f\t  %f\n", table[i].second.Z, table[i].second.N, table[i].second.l, evap.P().x0, evap.P().v3.x, evap.P().v3.y, evap.P().v3.z, table[i].second.j/2.0, table[i].second.E+random_E);
  //fclose(decay_file);

  products[0]=daughter; //updated
  products.push_back(evap); //new
}
