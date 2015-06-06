#include "separation_energy.hh"
#include "globals.hh"


double separation_E(Nucleus n)
{
  //get binding energy of neutrons and protons
  //by comparing mass of adjacent nuclei
  Nucleus nZm1;
  Nucleus nNm1;
  Nucleus N;
  Nucleus H1;
  nZm1.set_Z(n.Z()-1);
  nNm1.set_N(n.N()-1);
  nZm1.set_N(n.N());
  nNm1.set_Z(n.Z());

  N.set_Z(0);
  N.set_N(1);
  H1.set_Z(1);
  H1.set_N(0);

  double BEn=mass_model.mass(nNm1)+mass_model.mass(N)-mass_model.mass(n);
  double BEp=mass_model.mass(nZm1)+mass_model.mass(H1)-mass_model.mass(n);
  return BEn>BEp ? BEn : BEp;
}
