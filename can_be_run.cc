#include "globals.hh"

bool can_be_run(int Z, int N)
{
  //checks the masses needed by model to determine pairing energies.
  return (mass_model.mass_exists(Z,N) &&  mass_model.mass_exists(Z+1,N) &&  mass_model.mass_exists(Z-1,N) &&  mass_model.mass_exists(Z,N+1) &&  mass_model.mass_exists(Z,N-1));
  
}
