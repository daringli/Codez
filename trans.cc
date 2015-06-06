#include <cstdlib>     // atoi, abs
#include <iostream>
#include <vector>
#include <cstdio>
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
  std::string Rrange;
  std::string Lrange;
  std::vector<std::string> Evap;
  int Lstart=0;
  int Lend=5;
  double Rstart=0;
  double Rend=20;
  double dr=0.1;
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("energy,E", po::value<std::string>(&Rrange), "The range E1-E2 to evaluate the transmission probability for. Default: 0 to 20 MeV.") 
    ("dE,d", po::value<double>(&dr), "The discretization to use. Default: 0.1 MeV")
    ("L,l", po::value<std::string>(&Lrange), "The range of orbital angular momentum l1-l2 to evaluate the transmission coefficient for. Default: 0-5") 
    ("evap,e", po::value<std::vector<std::string> >(&Evap)->multitoken(), "The particles to evaporate by. Default: all. Options: p,n,d,t,He3,alpha. Additional particles can be invoked by using the flag again, i.e. -e He3 -e alpha.") 
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  if (vm["evap"].empty()) {
    //by default, we run over all default evaporations.
    Evap.push_back("p");
    Evap.push_back("n");
    Evap.push_back("d");
    Evap.push_back("t");
    Evap.push_back("He3");
    Evap.push_back("alpha");
  }

  if(vm.count("energy")){
    int sep=Rrange.find("-");
    if(sep==std::string::npos){
      Rstart=atof(Rrange.c_str());
      Rend=Rstart;
    }
    else{
      Rstart=atof(Rrange.substr(0,sep).c_str());
      Rend=atof(Rrange.substr(sep+1,std::string::npos).c_str());
    }
  }

  if(vm.count("L")){
    int sep=Lrange.find("-");
    if(sep==std::string::npos){
      Lstart=atof(Lrange.c_str());
      Lend=Lstart;
    }
    else{
      Lstart=atof(Lrange.substr(0,sep).c_str());
      Lend=atof(Lrange.substr(sep+1,std::string::npos).c_str());
    }
  }


  //initialize model for the transmission
  Codex_particle_model model;
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
    Nucleus f;
    Nucleus evap;
    //print info in comments
    printf("#Z=%i\tN=%i\n",Z,N);
    for(int i=0;i<Evap.size();i++){
      if(Evap[i]=="p"){
	f.set_N(N);
	f.set_Z(Z-1);
	evap.set_N(0);
	evap.set_Z(1);	
      }
      else if(Evap[i]=="n"){
	f.set_N(N-1);
	f.set_Z(Z);
	evap.set_N(1);
	evap.set_Z(0);
      }
      else if(Evap[i]=="d"){
	f.set_N(N-1);
	f.set_Z(Z-1);
	evap.set_N(1);
	evap.set_Z(1);
      }
      else if(Evap[i]=="t"){
	f.set_N(N-2);
	f.set_Z(Z-1);
	evap.set_N(2);
	evap.set_Z(1);
      }
      else if(Evap[i]=="He3"){
	f.set_N(N-1);
	f.set_Z(Z-2);
	evap.set_N(1);
	evap.set_Z(2);
      }
      else{
	f.set_N(N-2);
	f.set_Z(Z-2);
	evap.set_N(2);
	evap.set_Z(2);
      }

      printf("#Zevap=%i\tNevap=%i\n",evap.Z(),evap.N());

      printf("#E\t");

      for(int l=Lstart; l<=Lend; l=l+1){
	printf("l=%i\t",l);
      }
      printf("\n");
      //print potential to stdout
      for(float r=Rstart; r<=Rend; r=r+dr){
	printf("%.2f\t",r);
	for(int l=Lstart; l<=Lend; l=l+1){
	  Proximity_potential pot=model.potential_properties(n,evap);
	  double value=model.transmission(n,f,r,l,pot);
	  printf("%e\t",value);
	}
	printf("\n");
      }
    }
  }  
}
