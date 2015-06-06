#include <stdlib.h>     // atoi, abs
#include <vector>
#include <stdio.h>
#include <fstream>
#include <string>
//to parse program options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>



#include "globals.hh"
#include "nucleus.hh"
#include "codex_particle_model.hh"


int main(int argc, char *argv[])
{
  std::string Erange;
  std::string Jrange;
  double Jstart,Jend;
  double Estart=0.5;
  double Eend=20;
  double dE=0.1;
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("energy,E,e", po::value<std::string>(&Erange), "The energy range E1-E2 to evaluate the level density for. Default: 0.5 to 20 MeV.") 
    ("dE,D,d", po::value<double>(&dE), "The energy discretization to use. Default: 0.1 MeV")
    ("spin,J", po::value<std::string>(&Jrange), "The spin range J1-J2 to evaluate the level density for. Default: not used, see --use-spin") 
    ("use-spin", "If specified, uses the spin of the nucleus in the list. Overrides other spin options.") 
    ("use-energy", "If specified, uses the energy of the nucleus in the list. Overrides other energy options.") 
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  //parse input
  bool useE=0;
  bool useJ=1; //default 1
  if(vm.count("use-energy")){
    useE=1;
  }

  if(vm.count("energy") && useE==0){
    int sep=Erange.find("-");
    if(sep==std::string::npos){
      Estart=atof(Erange.c_str());
      Eend=Estart;
    }
    else{
      Estart=atof(Erange.substr(0,sep).c_str());
      Eend=atof(Erange.substr(sep+1,std::string::npos).c_str());
    }
  }

  if(vm.count("spin") && (vm.count("use-spin")==0)){
    useJ=0; //overrides default useJ=1
    int sep=Jrange.find("-");
    if(sep==std::string::npos){
      Jstart=atof(Jrange.c_str());
      Jend=Jstart;
    }
    else{
      Jstart=atof(Jrange.substr(0,sep).c_str());
      Jend=atof(Jrange.substr(sep+1,std::string::npos).c_str());
    }
  }


  //initialize model for level densities
  Codex_particle_model rhomodel;
  //get nuclei to run from stdin
  int Z,N;
  double E,J,px,py,pz;

  char line[500];
  char trimed_line[500];
  //scan stdin for data about the nuclei
  while(fgets(line, 500, stdin)){
    sscanf(line, "%s", trimed_line); //removes leading whitespace
    //ignore comment and event number lines
    if(line[0]=='*' || line[0]=='#'){
      continue;
    }  
    sscanf(line, "%d%d%lf%lf%lf%lf%lf", &Z,&N,&E,&J,&px,&py,&pz);
    Nucleus n;
    n.set_Z(Z);
    n.set_N(N);
    if(useE){
      Estart=E;
      Eend=E;
    }
    if(useJ){
      Jstart=J;
      Jend=J;
    }
    //print info in comments
    printf("#Z=%i\tN=%i\n",Z,N);
    printf("#E\t");
    for(float J=Jstart; J<=Jend; J=J+0.5){
      printf("J=%.1f\t",J);
    }
    printf("\n");
    //print level density to stdout
    for(float E=Estart; E<=Eend; E=E+dE){
      printf("%.2f\t",E);
      for(float J=Jstart; J<=Jend; J=J+0.5){
	n.set_E(E);
	n.set_J(J);
	double rho=rhomodel.rho(n,E);
	printf("%e\t",rho);
      }
      printf("\n");
    }
  }  
}
