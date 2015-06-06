#include "print_products.hh"
#include <cstdio>

void print_products(std::vector<std::vector<Nucleus> > ns, FILE* fout)
{
  //loop through all events
  for(int i=0; i<ns.size();i++){
    fprintf(fout,"*Event %i\n",i+1);
    //loop through the particles of the event
    for(int j=0; j<ns[i].size();j++){ 
      fprintf(fout,"%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",ns[i][j].Z(),ns[i][j].N(),ns[i][j].E(),ns[i][j].J(),ns[i][j].P().v3.x,ns[i][j].P().v3.y,ns[i][j].P().v3.z);
    }
  }
}
