#include <stdio.h> //printf
#include <string.h> //string
#include <fstream>
#include <iostream>
#include <sstream>

#include "globals.hh"
#include "mass_model.hh"

bool Mass_model::mass_exists(int Z, int N) const
{
  //returns true if the mass exists in the table.
  int A=Z+N;
  return masses.count(std::make_pair(Z,A));
 
}

double Mass_model::excess_mass(Nucleus n) const
{
  //returns the mass EXCESS of a given nuclei, i.e. it's deviation from A*(atomic mass unit). values in MeV
  int A=n.A();
  int Z=n.Z();
  if(A!=0){
    return masses.at(std::make_pair(Z,A));
  }
  else{
    return 0;
  }
}

double Mass_model::mass(Nucleus n) const
{
  //returns the mass EXCESS of a given nuclei, i.e. it's deviation from A*(atomic mass unit). values in MeV
  return excess_mass(n)+u*n.A();
}

void Mass_model::read_mass_file(const char*& MassDataFilename)
{
  int Z, Amax, Amin;
  double mass;
  int i=1;
  int j;
  char line[500];
  FILE * massFile;
  massFile = fopen(MassDataFilename,"r");
  //massFile = fopen("data/Masswbn.dat","r");

  //loop through every line in file
  while(fgets(line, 500, massFile))
    {
      //for odd numbered lines, get Z Amin, Amax
      if(i%2){
	sscanf(line, "%d%d%d", &Z,&Amin,&Amax);
      }
      else{
	//printf("Z: %i\n",Z);

	j=0;
	//for even numbered lines, get mass of nuclei by using stringstream magic.
	std::string linestring=std::string(line);
	std::stringstream linestream(linestring);
	while(linestream>>mass)
	  {
	    masses[std::make_pair(Z,Amin+j)] = mass;
	    //printf("%d mass in constructor:%f\n",j,mass);
	    //printf("%d mass in constructor:%f\n",j,masses[std::make_pair(Z,Amin+j)]);
	    j++;
	  }
      }
      i++;   
    }
  fclose(massFile);
}
