#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "globals.hh"

#include "decay_chain.hh"
#include "file_names.hh"


int main(int argc, char *argv[])
{

  if(argc!=2){
    printf("Error, incorrect number of arguments. Usage: %s input_file\n", argv[0]);
    return 1;
  }
  int num_mothers=10000;

  int lmax=10;
  int jmax=21;
  //get nuclei to run from file
  int Z,N;
  double E,J;

  FILE* input_file;
  input_file = fopen(argv[1],"r");
  char line[500];
  //scan file for lines on the form Z N E J
  while(fgets(line, 500, input_file))
    {
      sscanf(line, "%d%d%lf%lf", &Z,&N,&E,&J);
      Nucleus n;
      n.set_Z(Z);
      n.set_N(N);
      n.set_E(E);
      n.set_J(J);
      //open file and write decay spectrum to it
      char* filename=outfile_name(Z,N,E,J);
      FILE* fout=fopen(filename,"w");
      fprintf(fout,"#Zev\t Nev\t l\t |Jf\t  Ef\n");
      decay_chain(n,fout,lmax,jmax,num_mothers);
      fclose(fout);

      //run TALYS on same nuclei

      //create talys input file
      FILE* f_talys=fopen("talys_in","w");
      fprintf(f_talys,"projectile 0\nelement    %i\nmass       %i\nenergy     efile\noutmain n\noutspectra y",Z,Z+N);
      fclose(f_talys);
      //create talys energy file
      FILE *efile =fopen("efile","w");
      fprintf(efile,"2 0 1 %f\n%f 0.5 \n%f 0.5",E,E,E);
      fclose(efile);
      //generate talys command line, put it in bash file and run it
      //command appends output to previous decay spectrum.
      char system_string[40]="";
      strcat(system_string,"talys < talys_in >> ");
      char talys_filename[40]="talys";
      strcat(talys_filename,filename);
      strcat(system_string,talys_filename);

      FILE* run_talys=fopen("run_talys.sh","w");
      fprintf(f_talys,"#!/bin/sh\n%s",system_string);
      fclose(run_talys);
     
      system("chmod u+x ./run_talys.sh");
      printf("Running:\n%s\n",system_string);
      system("./run_talys.sh");
    }
  fclose(input_file);

}
