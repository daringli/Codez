#include <cmath> //pow
#include <boost/random/mersenne_twister.hpp> //randomness

#include "mass_model.hh"
#include <ctime>

boost::random::mt19937_64 gen(std::time(0));

//data with masses for nuclei
const char* mass_data="data/Masswbn.dat";
Mass_model mass_model(mass_data);

//physical constants
double e2=1.4399764; //elementary charge squared [MeV fm]
double hbar=197.0; //[MeV fm]
double hbar2=pow(hbar,2);
double u=931.5; //u in MeV

//energy discretization
double dE=0.1; //MeV

