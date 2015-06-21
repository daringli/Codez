#ifndef GLOBALS__HH
#define GLOBALS__HH

#include <boost/random/mersenne_twister.hpp> //randomness

#include "mass_model.hh"

//PRGN generator 
extern boost::random::mt19937_64 gen; 

//data with masses for nuclei
extern Mass_model mass_model;

//physical constants
extern double e2; //elementary charge squared [MeV fm]
extern double hbar; //in fm MeV
extern double hbar2;
extern double u; //u in MeV
extern double r0; //nuclear radius constant in fm


//energy discretization
extern double dE; //MeV

//mathematical constants
//no math here! M_PI used in code, 
//but if we decide to conform to standard, we may want to change that.

#endif
