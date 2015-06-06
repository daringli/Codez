#include <cstdlib>     // atoi, abs
#include <fstream>
#include <string>
#include <iostream>

//to parse program options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "can_be_run.hh"

//static int verbose_flag;

int main(int argc, char *argv[])
{
  double E,J,px,py,pz;
  int Zstart, Zend;
  int Nstart, Nend;
  int Astart, Aend;

  std::string Zrange,Nrange,Arange,Estring, Jstring;
  //to deal with mutual exclusive A,N,Z combinations
  bool use_A=0;
  bool use_N=0;
  bool use_Z=0;
 
  namespace po = boost::program_options; 
  po::options_description desc("Options"); 

  desc.add_options() 
    ("help,h", "Print help messages") 
    ("energy,E,e", po::value<std::string>(&Estring), "Energy to use.") 
    ("spin,J", po::value<std::string>(&Jstring), "Spin to use.") 
    ("Z,Z,z", po::value<std::string>(&Zrange), "Z range to use. Can be given Z1-Z2 or just Z. Z1>Z2. Z should be an integer") 
    ("N,N,n", po::value<std::string>(&Nrange), "N range to use. Can be given N1-N2 or just N. N1>N2. N should be an integer") 
    ("A,A,a", po::value<std::string>(&Arange), "A range to use. Can be given A1-A2 or just A. A1>A2. A should be an integer") 
    ("px,x", po::value<double>(&px), "Nucleus momentum in the x direction.") 
    ("py,y", po::value<double>(&py), "Nucleus momentum in the y direction.") 
    ("pz,z", po::value<double>(&pz), "Nucleus momentum in the z direction.") 
    ;
 
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);  

  //if help, print description.
  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  if ( vm.count("energy")==0) {  
    printf("No energy given! Use -E <energy> and try again.\n");
    return 1;
  }
  if ( vm.count("spin")==0) {  
    printf("No spin given! Use -J <spin> and try again.\n");
    return 1;
  }

  E=atof(Estring.c_str());
  if(isnan(E)){
    //parse option string
    printf("Not supported!\n");
  }
  //printf("E: %f\n",E);

  J=atof(Jstring.c_str());
  if(isnan(J)){
    //parse option string
    printf("Not supported!\n");
  }
  //printf("J: %f\n",J);
  
 
  if ( vm.count("Z")) { 
    int sep=Zrange.find("-");
    if(sep==std::string::npos){
      Zstart=atoi(Zrange.c_str());
      Zend=Zstart;
    }
    else{
      Zstart=atoi(Zrange.substr(0,sep).c_str());
      Zend=atoi(Zrange.substr(sep+1,std::string::npos).c_str());
    }
    //printf("Z in [%i,%i]\n",Zstart,Zend);
    use_Z=1;
  } 

if ( vm.count("N")) { 
    int sep=Nrange.find("-");
    if(sep==std::string::npos){
      Nstart=atoi(Nrange.c_str());
      Nend=Nstart;
    }
    else{
      Nstart=atoi(Nrange.substr(0,sep).c_str());
      Nend=atoi(Nrange.substr(sep+1,std::string::npos).c_str());
    }
    //printf("N in [%i,%i]\n",Nstart,Nend);
    use_N=1;
  } 

if ( vm.count("A")) { 
    int sep=Arange.find("-");
    if(sep==std::string::npos){
      Astart=atoi(Arange.c_str());
      Aend=Astart;
    }
    else{
      Astart=atoi(Arange.substr(0,sep).c_str());
      Aend=atoi(Arange.substr(sep+1,std::string::npos).c_str());
    }
    //printf("A in [%i,%i]\n",Astart,Aend);
    use_A=1;
  } 

//3 input Z,N,A
 if(use_A && use_N && use_Z){
   printf("Too many arguments given: A, N and Z.\n");
   return 2;
 }
 //2 input Z,N,A
 int first_start, second_start, first_end, second_end;
 if(use_A && use_N){
   first_start=Nstart;
   second_start=Astart;
   first_end=Nend;
   second_end=Aend;
 }
 if(use_A && use_Z){
   first_start=Zstart;
   second_start=Astart;
   first_end=Zend;
   second_end=Aend;  
 }
 if(use_N && use_Z){
   first_start=Zstart;
   second_start=Nstart;
   first_end=Zend;
   second_end=Nend;  
 }

  
 //We now generate a list of eligible nuclei given the input
 //first print a comment describing the output
 printf("#Z\tN\tE\tJ\tpx\tpy\tpz\n");

  
  for(int first=first_start;first<=first_end;first++){
    int i=1;
    for(int second=second_start;second<=second_end;second++){
      //the different cases:
      int Z,N;
      if(use_A && use_N){
	N=first;
	Z=second-first;
      }
      if(use_Z && use_A){
	Z=first;
	N=second-first;
      }
      if(use_Z && use_N){
	Z=first;
	N=second;
      }
      if(can_be_run(Z,N)){
	printf("*Event %i\n",i);
	printf("%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",Z,N,E,J,px,py,pz);
	i++;
      }
    }
  }
      
}
