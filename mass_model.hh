#ifndef MASS_MODEL__HH
#define MASS_MODEL__HH

#include <map>

#include "nucleus.hh"

class Mass_model
{
public:
  std::map<std::pair<int,int>, float> masses;


  Mass_model(const char*& MassDataFilename)
  {
    read_mass_file(MassDataFilename);
  }

  double excess_mass(Nucleus n) const;

  bool mass_exists(int Z, int N) const;

  double mass(Nucleus n) const;

  void read_mass_file(const char*& MassDataFilename);
};
#endif
