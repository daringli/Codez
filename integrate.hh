#ifndef INTEGRATE__HH
#define INTEGRATE__HH

#include "potential.hh"
#include <cstdlib>     //abs
double integrate(Potential const & pot, double E, double ri, double ro);


#endif
