#include "print_gunfile.hh"
#include "elements.hh"


void print_gunfile(std::vector<std::vector<Nucleus> > ns, FILE* fout)
{
  //loop through all events
  for(int i=0; i<ns.size();i++){
    fprintf(fout,"*** EVENT %i ***\n",i);
    fprintf(fout,"origin: {\n");
    //loop through the particles of the event
    for(int j=0; j<ns[i].size();j++){ 
      fprintf(fout,"  %s:",element_name(ns[i][j]));	
      fprintf(fout,"     pxyz=(%e, %e, %e)\n",ns[i][j].P().v3.x,ns[i][j].P().v3.y,ns[i][j].P().v3.z);
    }
    fprintf(fout,"}\n");
  }
}
