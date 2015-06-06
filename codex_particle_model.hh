#ifndef CODEX_PARTICLE_MODEL__HH
#define CODEX_PARTICLE_MODEL__HH

#include "model.hh"
#include "potential.hh"
#include <vector>
#include <exception>


struct no_minima_exception: public std::exception
{
  virtual const char* what() const throw()
  {
    return "No minima in potential found";
  }
};

struct xFx{
  //stores x and F(x).
  double f;
  double x;
};

class Proximity_potential: public Potential
{
private:
  int m_l;
  double m_rin;
  double m_coul1;
  double m_coul2;
  double m_rc;
  double m_b; //slope of level density at fermi level.
  double m_vnConst;
  double m_mu;
  double m_factor;

public:
  Proximity_potential(double rinSet, double coulSet, double rcSet, double BSet,double vnConstSet, double muSet);

  double value(double r) const;

  double first_derivative(double r) const;

  double second_derivative(double r) const;

  void print_to_file(char* filename);

  void print_derivative_to_file(char* filename);

  double rc()const{return m_rc;}
  double rin()const{return m_rin;}
  double vn_const()const{return m_vnConst;}
  double mu()const{return m_mu;}
  int l()const{return m_l;}

  //may want to add more getters... and setters!
  void set_l(int L){m_l=L;}

};

class Codex_particle_model : public Model
{
public:
  //Codex_particle_model(const char*& MassDataFilename)
  //: Model(MassDataFilename){}; //inherits constructor.

  void Rsum(Nucleus n, int jmax,int lmax,std::vector<std::pair<double,Decay> > & Rsum_decay,  int & counter) const;

  Proximity_potential potential_properties(Nucleus mother, Nucleus evaporation) const;

  std::vector<xFx> potential_min_max(Proximity_potential const & pot, double x, int numChanges) const;

  double tunneling(double E, Proximity_potential pot, xFx max, xFx min) const;

  double transmission(Nucleus initial, Nucleus final, double E, int l, Proximity_potential pot) const;

  double spin(Nucleus n) const;
  
  virtual ~Codex_particle_model(){};
};



#endif
