#ifndef MODEL__HH
#define MODEL__HH

#include <vector>
#include "nucleus.hh"
#include "decay.hh"

struct Primary_inertia
{
  double x;
  double y;
  double z;
};


class Model
{
public:

  virtual void Rsum(Nucleus n, int jmax,int lmax,std::vector<std::pair<double,Decay> > & Rsum_decay,  int & counter) const =0;
    
  double macroscopic_mass(Nucleus n) const;
    
  double gs_deformation_energy(Nucleus n) const;

  double microscopic_mass(Nucleus n) const;

  double gs_beta(Nucleus n) const;

  Primary_inertia primary_inertia(Nucleus n) const;

  double yrast(Nucleus n) const;

  double level_density_parameter(Nucleus n) const;

  double pairing(Nucleus n) const;

  double critical_energy(Nucleus n) const;

  double shell_damping(Nucleus n, double Eeff) const;

  double pair_damping(Nucleus n,double Eeff) const;

  double rho(Nucleus n, double E) const;

  double intrinsic_energy(Nucleus n) const;

  virtual ~Model(){};
};


#endif
