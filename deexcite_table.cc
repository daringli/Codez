#include <cstdio>

#include "deexcite_table.hh"

std::vector<std::pair<double,Decay> > deexcite_table(Nucleus mother, const std::vector<Model*> &models,int jmax, int lmax)
{
  std::vector<std::pair<double,Decay> > Rsum_Decay(1);
  Rsum_Decay[0].first=0;
  Decay decay;
  //0 is empty since we add sum in Rsum_Decay.first, which we need to zero initialize, 
  //but we can't just push new elements in std::vector if 0 isn't fully filled.
  Rsum_Decay[0].second=decay; 
  int counter=0;
  int num_models=models.size();
  for(int i; i<num_models;i++){
    models[i]->Rsum(mother, jmax,lmax,Rsum_Decay, counter);
  }
  //cannot use resize since no default constructor. erase instead
  Rsum_Decay.erase(Rsum_Decay.begin()+counter+1,Rsum_Decay.begin()+Rsum_Decay.size());
  return Rsum_Decay;
}
