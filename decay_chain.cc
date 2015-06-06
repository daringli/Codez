#include <cstdlib>     // atoi, abs
#include <vector>
#include <cstdio>
#include <fstream>
#include <ctime>

#include <boost/random/mersenne_twister.hpp> //randomness

#include "nucleus.hh"
#include "model.hh"
#include "codex_particle_model.hh"
#include "codex_gamma_model.hh"
#include "deexcite_table.hh"
#include "deexcite.hh"
#include "quasi_elastic.hh"
#include "print_gunfile.hh"
#include "globals.hh"

void decay_chain(std::vector<Nucleus> &nuclei,const std::vector<Model*> & models, int lmax, int jmax, int num_decays, double  E_threshold)
{
  //do decay steps until energy is low enough.

 
  int decay_counter=0;
  //do decay steps until energy is low enough.
  while(nuclei[0].E()>E_threshold && decay_counter<num_decays){
    //first generate table with transition probabilities
    //ONLY FOR NUCLEI[0] at the moment
    std::vector<std::pair<double,Decay> > decay_table=deexcite_table(nuclei[0], models,jmax,lmax);

    //deexcite with table
    deexcite(nuclei,decay_table);
    decay_counter++;
  }
}
