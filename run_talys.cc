#include <cstring>
#include <cstdlib>
#include <string>
#include <cstdio>
#include "nucleus.hh"

void run_talys(Nucleus n, std::string infile_name, std::string outfile_name)
{
  //generates a TALYS input-file, and an energy file to specify an
  //initial population with a single energy¨
  //NOTE: requires TALYS 1.73, "released" 07-05-2015.
  int Z=n.Z();
  int A=n.A();
  double E=n.E();

  FILE* f_intalys=fopen(infile_name.c_str(),"w");
  fprintf(f_intalys,"projectile 0\nelement    %i\nmass       %i\nenergy     efile\noutmain n\noutspectra y \noutbinspectra y\npopmev y\nmaxz 0\nmaxn 0",Z,A);
  fclose(f_intalys);
  //create talys energy file
  FILE *efile =fopen("efile","w");
  fprintf(efile,"4 0 1 \n0.0 0.0 \n%f 0.0\n%f 1.0\n%f 0.0",E-0.5,E,E+0.5);
  fclose(efile);
  //generate talys command line, put it in bash file and run it
  //command appends output to previous decay spectrum.
  char system_string[40]="";
  strcat(system_string,"talys1-73 < ");
  strcat(system_string,infile_name.c_str());
  strcat(system_string," > ");
  strcat(system_string,outfile_name.c_str());
 
  FILE* run_talys=fopen("run_talys.sh","w");
  fprintf(run_talys,"#!/bin/sh\n%s",system_string);
  fclose(run_talys);
     
  system("chmod u+x ./run_talys.sh");
  //printf("# %s\n",system_string);
  system("./run_talys.sh");
}
