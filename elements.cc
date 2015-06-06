#include <cstring>
#include <cstdio>
#include "elements.hh"

char* element_name(Nucleus n)
{
  int Z=n.Z();
  int N=n.N();
  if(Z==0){
    if(N==1)
      return "n";
    if(N==0)
      return "gamma";      
  }

  //Such beauty...
  const char *elements[] =
    {
      "H",                                                                                 "He",
      "Li","Be",                                                  "B", "C", "N", "O", "F", "Ne",
      "Na","Mg",                                                  "Al","Si","P", "S", "Cl","Ar",
      "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
      "Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I", "Xe",
      "Cs","Ba","La",
      /**/           "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
      /**/           "Hf","Ta","W", "Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
      "Fr","Ra","Ac",
      /**/           "Th","Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
      /**/           "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Ed","Fl","Ef","Lv","Eh","Ei",
    };

  static char ret[10];
  sprintf(ret,"%s%i",elements[Z-1],n.A());
  return ret;
}
