#ifndef DEEXCITE__HH
#define DEEXCITE__HH

#include <vector>
#include "nucleus.hh"
#include "decay.hh"
#include <vector>

void deexcite(std::vector<Nucleus> &products, const std::vector<std::pair<double,Decay> > &table);
		
#endif
