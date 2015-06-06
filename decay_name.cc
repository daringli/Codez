#include "decay_name.hh"

std::string decay_name(Decay d){
  if(d.Z==0 && d.N==0)
    return "gamma";
  if(d.Z==1 && d.N==0)
    return "p";
  if(d.Z==0 && d.N==1)
    return "n";
  if(d.Z==1 && d.N==1)
    return "d";
  if(d.Z==1 && d.N==2)
    return "t";
  if(d.Z==2 && d.N==1)
    return "He3";
  if(d.Z==2 && d.N==2)
    return "alpha";
  else return "ERROR"; //this should never happen
}
