#include <cstdlib>     // atoi, abs
#include <fstream>
#include <string>
#include <iostream>

//to parse program options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "globals.hh"
#include "nucleus.hh"

#include "decay_chain.hh"
#include "lorentz_boost.hh"
#include "vector34.hh"
#include "quasi_elastic.hh"
#include "codex_particle_model.hh"
#include "codex_gamma_model.hh"
#include "deexcite_table.hh"
#include "deexcite.hh"
#include "decay_chain.hh"
#include "print_gunfile.hh"

//static int verbose_flag;

int main(int argc, char *argv[])
{  
  int proj_Z;
  int proj_N;
  double proj_T;

  int target_Z;
  int target_N;

  int cluster_Z;
  int cluster_N;

  int num_events=10; //conservative default!

  std::string Estring, Jstring;

  namespace po = boost::program_options; 
  po::options_description desc("Options"); 

  desc.add_options() 
    ("help,h", "Print help messages") 
    ("ex-energy,E", po::value<std::string>(&Estring), "Excitation energy of prefragment. [MeV]") 
    ("spin,J", po::value<std::string>(&Jstring), "Spin of prefragment.") 
    ("kinetic,T", po::value<double>(&proj_T), "Kinetic energy of the projectile  [MeV].") 
    ("events,N", po::value<int>(&num_events), "Number of events to generate.") 
    ("target-Z", po::value<int>(&target_Z), "Z of the target.") 
    ("target-N", po::value<int>(&target_N), "N of the target.") 
    ("proj-Z", po::value<int>(&proj_Z), "Z of the projectile.") 
    ("proj-N", po::value<int>(&proj_N), "N of the projectile.") 
    ("cluster-Z", po::value<int>(&cluster_Z), "Z of the cluster.") 
    ("cluster-N", po::value<int>(&cluster_N), "N of the cluster.") 
    ;


  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);  

  //if help, print description and quit.
  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  //if projectile kinetic energy is not specified
  if ( vm.count("kinetic")==0) {  
    printf("No beam kinetic energy given! Use -T <energy> and try again.\n");
    return 1;
    //possibly use a default model
  }

  //if prefragment is not specified
  if ( vm.count("ex-energy")==0) {  
    printf("No energy given! Use -E <energy> and try again.\n");
    return 1;
    //possibly use a default model
  }
  if ( vm.count("spin")==0) {  
    printf("No spin given! Use -J <spin> and try again.\n");
    return 1;
    //possibly use a default model
  }

  //if target, cluster or projectile not given
  if ( (vm.count("target-Z")==0) && (vm.count("target-N")==0)) {  
    printf("No target given! Use --target-Z <Z> --target-N <N> and try again.\n");
    return 1;
    //possibly use a default model
  }
  if ( (vm.count("proj-Z")==0) && (vm.count("proj-N")==0)) {  
    printf("No projectile given! Use --proj-Z <Z> --proj-N <N> and try again.\n");
    return 1;
    //possibly use a default model
  }
  if ( (vm.count("cluster-Z")==0) && (vm.count("cluster-N")==0)) {  
    printf("No cluster given! Use --cluster-Z <Z> --cluster-N <N> and try again.\n");
    return 1;
    //possibly use a default model
  }

  //at this point, we know that the input is somewhat sane, so we parse it!

  double prefragment_excitation_energy=atof(Estring.c_str());
  if(isnan(prefragment_excitation_energy)){
    //parse option string
    printf("Not supported!\n");
    return 1;
  }
  //printf("E: %f\n",E);

  double prefragment_J=atof(Jstring.c_str());
  if(isnan(prefragment_J)){
    //parse option string
    printf("Not supported!\n");
    return 1;
  }
  //printf("J: %f\n",J);
  
  //declare nucleus objects

  //declare projectile object
  Nucleus projectile;
  projectile.set_Z(proj_Z);
  projectile.set_N(proj_N);
 double proj_mass=mass_model.mass(projectile);
  Vector4 proj_P;
  proj_P.x0=proj_T+proj_mass;
  proj_P.v3.x=0;
  proj_P.v3.y=0;
  proj_P.v3.z=sqrt(pow(proj_P.x0,2)-pow(proj_mass,2));
  projectile.set_P(proj_P);
 
  //declare target object
  Nucleus target;
  target.set_Z(target_Z);
  target.set_N(target_N);

  //declare cluster object
  Nucleus cluster;
  cluster.set_Z(cluster_Z);
  cluster.set_N(cluster_N);

  //vector of reaction products
  std::vector<Nucleus> products;

  //generate+print results of several quasi-elastic scatterings with prefragment with set J, E*
  printf("#Z\tN\tE\tJ\tpx\tpy\tpz\n");
  for(int i=1; i<=num_events;i++){
    printf("*Event %i\n",i);
    std::vector<Nucleus> products=quasi_elastic(projectile,target,cluster,prefragment_excitation_energy,prefragment_J);

    //print info about this event.
    int num_prods=products.size();
    //printf("num_prods: %i\n",num_prods);
    for(int j=0;j<num_prods;j++){
      //printf("j:%i\n",j);
      printf("%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",products[j].Z(),products[j].N(),products[j].E(),products[j].J(),products[j].P().v3.x,products[j].P().v3.y,products[j].P().v3.z);
    }    
  }      
}
