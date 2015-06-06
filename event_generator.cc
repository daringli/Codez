#include <cstdlib>     // atoi, abs
#include <cmath>
#include <ctime>
#include <fstream>

#include <boost/random/mersenne_twister.hpp> //randomness
#include <boost/random/uniform_real_distribution.hpp> //uniform randomness!
#include <boost/random/normal_distribution.hpp>

#include "globals.hh"

#include "decay_chain.hh"
#include "lorentz_boost.hh"
#include "vector34.hh"
#include "quasi_elastic.hh"
#include "codex_particle_model.hh"
#include "codex_gamma_model.hh"
#include "deexcite_table.hh"
#include "deexcite.hh"
#include "decay_chain.hh"
#include "print_gunfile.hh"


int main(int argc, char *argv[])
{

  if(argc!=2){
    printf("Error: incorrect number of input arguments\n Usage: %s gunfile-name\n",argv[0]);
    return 1;
  }

  //calculation for the collision process
  //---------"user input" that decides what gets knocked out-----------
  //projectile
  int projectile_Z=9;
  int projectile_N=8;
  double projectile_T=500*(projectile_Z+projectile_N); //MeV

  //target
  int target_Z=1;
  int target_N=0;

  //cluster, the nucleons in the projectile that participates in the collision
  //pn scattering
  int cluster_Z=0;
  int cluster_N=1; 

 int num_events=100; //number of decay process to initiate for this particular prefragment.
 
  double prefragment_excitation_energy=20;
  double prefragment_J=0; //need model for J?

  //------Declare nuclei objects-------------

  //declare projectile object
  Nucleus projectile;
  projectile.set_Z(projectile_Z);
  projectile.set_N(projectile_N);
  double projectile_mass=mass_model.mass(projectile);
  Vector4 projectile_P;
  projectile_P.x0=projectile_T+projectile_mass;
  projectile_P.v3.x=0;
  projectile_P.v3.y=0;
  projectile_P.v3.z=sqrt(pow(projectile_P.x0,2)-pow(projectile_mass,2));
  projectile.set_P(projectile_P);
    
  //declare target object
  Nucleus target;
  target.set_Z(target_Z);
  target.set_N(target_N);

  //declare cluster object
  Nucleus cluster;
  cluster.set_Z(cluster_Z);
  cluster.set_N(cluster_N);
 
  //vector of vector of reaction products
  //right now only products[i][0] is assumed to be excited
  //the other elements are all the ejectiles from the event.
  std::vector<std::vector<Nucleus> > products(num_events);

  for(int i=0; i<num_events;i++){
    //generate product and prefragment with set J, E* by quasi-elastic sctring.
    //identical prefragments and products, but different momenta each event.
    products[i]=quasi_elastic(projectile,target,cluster,prefragment_excitation_energy,prefragment_J);
    //printf("qs i: %i\n",i);
  }
  //initialize model objects
  Codex_particle_model p_model;
  Codex_gamma_model g_model;
  std::vector<Model*> models(2);
  models[0]=&g_model;
  models[1]=&p_model;

  int lmax=10;
  int jmax=21; //convention, Jf=jf/2. Note that small js are integers. 
  
  //create table with transition probabilities from the initial prefragment
  std::vector<std::pair<double,Decay> > initial_decay_table=deexcite_table(products[0][0], models,jmax,lmax);

  for(int i=0; i<num_events;i++){
    printf("deexcite i: %i\n",i);

    //deexcite the initial prefragments with the same transition table. 
    //Update products with the next step in chain.
    deexcite(products[i],initial_decay_table);

    //we now have different states for each product, and thus need to generate 
    //a table for each event, which the decay_chain function does repeatedly
    //until the prefragment is deexcited.
    decay_chain(products[i], models, lmax, jmax);
  }

  //print all the products to a gunfile
  FILE* gunfile=fopen(argv[1],"w");
  print_gunfile(products,gunfile);
  fclose(gunfile);

  return 0;
}
