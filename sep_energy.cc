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


int main(int argc, char *argv[])
{
  std::string Rrange;
  std::string Lrange;
  std::vector<std::string> Evap;
  int Lstart=0;
  int Lend=5;
  double Rstart=0.5;
  double Rend=20;
  double dr=0.1;
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("part,p", po::value<std::vector<std::string> >(&Evap)->multitoken(), "The particles to calculate the separation energy for. Default: all. Options: p,n,d,t,He3,alpha. Additional particles can be invoked by using the flag again, i.e. -p He3 -p alpha.") 
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  if (vm["part"].empty()) {
    //by default, we run over all default evaporations.
    Evap.push_back("p");
    Evap.push_back("n");
    Evap.push_back("d");
    Evap.push_back("t");
    Evap.push_back("He3");
    Evap.push_back("alpha");
  }

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
    int Nev,Zev;

    printf("#Z=%i\tN=%i\n",Z,N);
    for(int i=0;i<Evap.size();i++){
      if(Evap[i]=="p"){
	Nev=0;
	Zev=1;
      }
      else if(Evap[i]=="n"){
	Nev=1;
	Zev=0;
      }
      else if(Evap[i]=="d"){
	Nev=1;
	Zev=1;
      }
      else if(Evap[i]=="t"){
	Nev=2;
	Zev=1;
      }
      else if(Evap[i]=="He3"){
	Nev=1;
	Zev=2;
      }
      else{
	Nev=2;
	Zev=2;
      }

      printf("#Zevap=%i\tNevap=%i\t",Zev,Nev);
      Nucleus ini;
      Nucleus fin;
      Nucleus ev;
      ini.set_N(N);
      ini.set_Z(Z);
      ev.set_N(Nev);
      ev.set_Z(Zev);
      fin.set_N(N-Nev);
      fin.set_Z(Z-Zev);

      double se=mass_model.excess_mass(ev)+mass_model.excess_mass(fin)-mass_model.excess_mass(ini);
      printf("%f\n",se);
      }
    }
}
