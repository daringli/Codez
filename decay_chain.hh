#ifndef DECAY_CHAIN__HH
#define DECAY_CHAIN__HH

#include "nucleus.hh"
#include "model.hh"

void decay_chain(std::vector<Nucleus> &nuclei,const std::vector<Model*> & models, int lmax, int jmax, int num_decays, double  E_threshold);

#endif
