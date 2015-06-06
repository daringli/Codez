#ifndef DEEXCITE_TABLE__HH
#define DEEXCITE_TABLE__HH

#include <vector>

#include "decay.hh"
#include "nucleus.hh"
#include "model.hh"


std::vector<std::pair<double,Decay> > deexcite_table(Nucleus mother, const std::vector<Model*> &models, int jmax,int lmax);

#endif
