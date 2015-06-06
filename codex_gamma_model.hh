#ifndef CODEX_PARTICLE_MODEL__H
#define CODEX_PARTICLE_MODEL__H

#include "model.hh"
#include <vector>

class Codex_gamma_model :public Model
{
public:
  //Codex_gamma_model(const char*& MassDataFilename)
  //: Model(MassDataFilename){}; //inherits constructor.

  void Rsum(Nucleus n, int jmax,int lmax,std::vector<std::pair<double,Decay> > & Rsum_decay, int & counter) const;

  double transmission(Nucleus initial, double E, int l) const;

};



#endif
