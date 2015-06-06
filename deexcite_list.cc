#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
//to parse program options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "globals.hh" //contains dE!!

#include "decay_chain.hh"
#include "lorentz_boost.hh"
#include "vector34.hh"
#include "quasi_elastic.hh"
#include "deexcite_table.hh"
#include "deexcite.hh"
#include "decay_chain.hh"
#include "codex_particle_model.hh"
#include "codex_gamma_model.hh"
#include "model.hh"
#include "print_products.hh"



int main(int argc, char *argv[])
{
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  double Jmax=10;
  int num_decays=99999999; //this is a huge default value.
  double E_threshold=5*dE; //stop when E* drops below this.
  int lmax=10;
  int jmax=21;
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("decays,d", po::value<int>(&num_decays), "Number of decays to preform. Default: unlimited") 
    ("threshold-E,e", po::value<double>(&E_threshold), "Excitation energy at which to stop decay. Default: 0.5MeV")
    ("dE", po::value<double>(&dE), "Energy discretization spacing. Default: 0.1MeV")
    ("max-J,J", po::value<double>(&Jmax), "The maximum spin to be considered for the nuclei.") 
    ("max-l,l", po::value<int>(&lmax), "The maximum orbital angular momentum to be considered for the emitted particle.")
    ;


 
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  if( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  if(vm.count("max-J")){
    jmax=2*Jmax;
  }

  //initialize models to use for decay
  Codex_particle_model pmodel;
  Codex_gamma_model gmodel;
  std::vector<Model*> models(2);
  models[0]=&gmodel;
  models[1]=&pmodel;

  //initilize event counter
  int event_counter=1;

  //get nuclei to run from stdin
  int Z,N;
  double E,J,px,py,pz;
  std::vector<std::vector<Nucleus> > all_products;
  std::vector<Nucleus> event_products;


  char line[500];
  char trimed_line[500];
  //scan stdin for data about the nuclei
  while(fgets(line, 500, stdin)){
    sscanf(line, "%s", trimed_line); //removes leading whitespace
    //ignore comments '#' and go to next event at '*'
    if(trimed_line[0]=='#'){
      continue;
    }  
    else if(trimed_line[0]=='*'){
      if(event_counter>1){
	//decay everything that was part of the previous event
	decay_chain(event_products, models, lmax, jmax, num_decays,E_threshold);
	all_products.push_back(event_products);
	event_products.clear();
      }
      event_counter++; //currently ignores the actual number that is printed
      continue;
    }  
    sscanf(line, "%d%d%lf%lf%lf%lf%lf", &Z,&N,&E,&J,&px,&py,&pz);
    Nucleus n;
    n.set_Z(Z);
    n.set_N(N);
    n.set_E(E);
    n.set_J(J);
    if(J>jmax/2.0){
      puts("Nuclei with J>Jmax! Aborting. Use -J to set a different Jmax");
      return 1;
    }
    //calculate the 0th component of the nucleus 4-momentum.
    Vector4 P;
    P.v3.x=px;
    P.v3.y=py;
    P.v3.z=pz;
    P.x0=sqrt(E*E+mass_model.mass(n)*mass_model.mass(n)+P.v3*P.v3);
    n.set_P(P);
    //add this nuclei to the current event.
    event_products.push_back(n);
  }
  //last event is added manually, since while(fgets) loop ended before '*'
  decay_chain(event_products, models, lmax, jmax, num_decays,E_threshold);
  all_products.push_back(event_products);

  //print decayed nuclei to stdout
  printf("#Z\tN\tE\tJ\tpx\tpy\tpz\n");
  print_products(all_products,stdout);
}

