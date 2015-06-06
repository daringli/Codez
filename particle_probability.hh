#ifndef PARTICLE_PROBABILITY__HH
#define PARTICLE_PROBABILITY__HH

#include <vector>
#include <string>
#include "decay.hh"
#include "nucleus.hh"

struct Particle_probability{
  std::string name;
  double prob;

  bool operator==(const Particle_probability& a) const
  {
    return (prob == a.prob);
  }

   bool operator<(const Particle_probability& a) const
  {
    return (prob < a.prob);
  }

   bool operator>(const Particle_probability& a) const
  {
    return (prob > a.prob);
  }

   bool operator<=(const Particle_probability& a) const
  {
    return (prob == a.prob) || (prob < a.prob);
  }

   bool operator>=(const Particle_probability& a) const
  {
    return (prob == a.prob) || (prob > a.prob);
  }

};

std::vector<Particle_probability> particle_width(std::vector<std::pair<double,Decay> > &table);

void print_table(std::vector<std::pair<double,Decay> > &table, Nucleus n);

#endif
