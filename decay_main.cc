#include <cstdlib>     // atoi, abs

#include "globals.hh"
#include "decay_chain.hh"
#include "nucleus.hh"

int main(int argc, char *argv[])
{

 if(argc!=6){
    printf("Error, incorrect number of arguments. Usage: %s Z N E* J outfile\n", argv[0]);
    return 1;
  }
  int Z=atoi(argv[1]);
  int N=atoi(argv[2]);
  double E=atof(argv[3]); //excitation energy [MeV]. 
  double J=atof(argv[4]);
  FILE* fout=fopen(argv[5],"w");

  Nucleus n;
  n.set_Z(Z);
  n.set_N(N);
  n.set_E(E);
  n.set_J(J);
  int num_mothers=1;


  int lmax=10;
  int jmax=21; //convention, Jf=jf/2. Note that small js are integers.
  fprintf(fout,"#Zev\t Nev\t l\t |Jf\t  Ef\n");
  decay_chain(n,fout,lmax,jmax,num_mothers);
  printf("Closing file %s\n----------\n",argv[5]);
  fclose(fout);
}
