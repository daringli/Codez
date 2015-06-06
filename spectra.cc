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


#include "globals.hh"
#include "vector34.hh"
#include "deexcite_table.hh"
#include "decay.hh"
#include "model.hh"
#include "codex_particle_model.hh"
#include "codex_gamma_model.hh"
#include "particle_probability.hh" 
#include "run_talys.hh"
#include "file_names.hh"
#include "talys_probability.hh"


int main(int argc, char *argv[])
{
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  std::string details;
  double Jmax;
  int lmax=10;
  int jmax=21;
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("details,d", po::value<std::string>(&details), "Level of detail of the output spectra. Possible options: full, particles, single; for all the details, the particle spectra, and just the most likely particle, respectively.") 
    ("talys,t", "Use Talys to generate the spectra") 
    ("max-J,J", po::value<double>(&Jmax), "The maximum spin to be considered for the nuclei.") 
    ("max-l,l", po::value<int>(&lmax), "The maximum orbital angular momentum to be considered for the emitted particle.") 
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  //if help, print description.
  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  bool to_run_talys=0;
  if(vm.count("talys")){
    to_run_talys=1;
  }

  if(vm.count("max-J")){
    jmax=2*Jmax;
  }
  //printf("jmax: %i\n",jmax);
  //printf("Jmax: %f\n",Jmax);
  //printf("lmax: %i\n",lmax);

  //initialize models to use for decay
  Codex_particle_model pmodel;
  Codex_gamma_model gmodel;
  std::vector<Model*> models(2);
  models[0]=&gmodel;
  models[1]=&pmodel;
  
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

    //printf("Z, N, E, J: %i, %i, %f, %f\n",Z,N,E,J); 
    //generate decay probabilities
    std::vector<std::pair<double,Decay> > table=deexcite_table(n, models,jmax,lmax);
    std::vector<Particle_probability> particle_widths;

    //sum up contributions of the particle decay channels in case the full spectrum isn't wanted
    if(details=="particles" || details=="single"){
      if(to_run_talys){
	std::string talys_infile="talys_in";
	std::string talys_outfile=outfile_name(Z,N,E,J,".out");
	run_talys(n,talys_infile,talys_outfile);
	particle_widths=talys_width(talys_outfile);
      }
      else{
	particle_widths=particle_width(table);
      }

      //greatest width last. First would be better...
      std::sort (particle_widths.begin(), particle_widths.end()); 

      //if we just want the most likely particle to decay by
      if(details=="single"){
	printf("%i\t%i\t%f\t%f\t%s\t%f\n",Z,N,n.E(),n.J(),particle_widths[particle_widths.size()-1].name.c_str(),particle_widths[particle_widths.size()-1].prob);
      }
      else if(details=="particles"){
	//if we want to see how likely every particle is
	for(int i=1;i<=particle_widths.size();i++){
	  printf("%i\t%i\t%f\t%f",Z,N,n.E(),n.J());
	  printf("\t%s\t%f\n",particle_widths[particle_widths.size()-i].name.c_str(),particle_widths[particle_widths.size()-i].prob);
	}
      }
    }
    else if(details=="full"){
      //puts("full!");
      print_table(table,n);
    }
  }
}
