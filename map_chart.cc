#include <cstdlib>     // atoi, abs
#include <cmath>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <string>


#include <boost/random/mersenne_twister.hpp> //randomness
#include <boost/random/uniform_real_distribution.hpp> //uniform randomness!
#include <boost/random/normal_distribution.hpp>

#include "globals.hh"
#include "nucleus.hh"

#include "decay_chain.hh"
#include "lorentz_boost.hh"
#include "vector34.hh"
#include "model.hh"
#include "codex_particle_model.hh"
#include "codex_gamma_model.hh"
#include "deexcite_table.hh"
#include "separation_energy.hh" 
#include "particle_probability.hh" 
#include "run_talys.hh"
#include "file_names.hh"
#include "talys_probability.hh"



int main(int argc, char *argv[])
{

    bool energy_given=0;
    bool to_run_talys=0;

  if(argc<2){
    printf("Error: incorrect number of input arguments\n Usage: %s outfile-name [Energy] [Parsed talys outfile]\n",argv[0]);
    return 1;
  }
  if(argc>2){
    energy_given=1;
    if(argc>3){
      to_run_talys=1;
    }
  }

  double J=0;

  //start of sweep of chart
  int Zlow=10;
  int Nlow=10; 
  //stop at
  int Zhigh=30;
  int Nhigh=30;
  //int Zhigh=30;
  //int Nhigh=30;

  int lmax=10;
  int jmax=21; //convention, Jf=jf/2. Note that small js are integers. 

  Codex_particle_model pmodel;
  Codex_gamma_model gmodel;
  std::vector<Model*> models(2);
  models[0]=&gmodel;
  models[1]=&pmodel;
  
  //crawl over the chart of the nuclides!
  FILE* fout=fopen(argv[1],"w");
  FILE* fouttalys;
  if(to_run_talys){
    fouttalys=fopen(argv[3],"w");
  }
  bool hit=0;
  for(int N=Nlow;N<=Nhigh;N++){
    for(int Z=Zlow;Z<=Zhigh;Z++){
      Nucleus n;
      double E;
      n.set_Z(Z);
      n.set_N(N);
      n.set_J(J);

      //determine separation energy and decay probabilities for this nucleus
      //will throw exception if mass not in table.
      std::vector<std::pair<double,Decay> > table;
      try 
	{
	  //need E*>>0 the model assumptions break down. 
	  //Want to be a bit above threshold for n and p, so they can compete.
	  if(energy_given){
	    E=atof(argv[2]);
	  }
	  else{
	    E=separation_E(n) >9 ? 1.1*separation_E(n) : 10.0; 
	  }
	  n.set_E(E);

	  printf("Z, N, E: %i, %i, %f\n",Z,N,E);

	  table=deexcite_table(n, models,jmax,lmax);
	} 
      catch (const std::out_of_range& oor)
	{
	  if(hit==0){
	    printf("next N!\n");
	    hit=1;
	    break;
	  }
	  else{
	    printf("try next Z!\n");
	    Zlow++;
	    continue;
	  }
	}
      hit=0;
      //generate partial width of each decay channel 
      std::vector<Particle_probability> particle_widths=particle_width(table);
      //greatest width last. First would be better...
      std::sort (particle_widths.begin(), particle_widths.end());
      printf("   by: %s with %f\%\n",particle_widths[particle_widths.size()-1].name.c_str(),particle_widths[particle_widths.size()-1].prob*100);
      fprintf(fout,"%i\t%i\t%f\t%s\t%f\n",Z,N,n.E(),particle_widths[particle_widths.size()-1].name.c_str(),particle_widths[particle_widths.size()-1].prob);

      if(to_run_talys){
	//run TALYS-1.73 for the same nuclei
	std::string talys_infile="talys_in";
	std::string talys_outfile=outfile_name(Z,N,E,J,".out");
	run_talys(n,talys_infile,talys_outfile);
	//parse TALYS output for particle_probabilities
	std::vector<Particle_probability> talys_widths=talys_width(talys_outfile);
	std::sort (talys_widths.begin(), talys_widths.end());
	fprintf(fouttalys,"%i\t%i\t%f\t%s\t%f\n",Z,N,n.E(),talys_widths[talys_widths.size()-1].name.c_str(),talys_widths[talys_widths.size()-1].prob);
      }
    }
    hit=1;
  }
  printf("Closing file %s\n----------\n",argv[1]);
  fclose(fout);
  if(to_run_talys){
    printf("Closing file %s\n----------\n",argv[3]);
    fclose(fouttalys);
  }
}
