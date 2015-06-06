#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm> //for find
#include <iostream> //for cout

//to parse program options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

int main(int argc, char *argv[])
{
  std::string out_filename;
  bool use_all=1; //by default, we use all the events on the list

  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("filename,f",  po::value<std::string> (&out_filename),"Filename to write to. ") 
    ;

 
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  //if help, print description.
  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

  if(vm.count("filename")==0){
    puts("Need a filename to write to!");
    return 1;
  }

  //prepare rootfile
  TFile f(out_filename.c_str(),"recreate");


  //get nuclei to run from stdin. 
  int Z,N,l;
  float deltaE,J,p;

  //prepare the tree and its branches
  TTree *tree = new TTree("spectra_out","Full spectra output");
  TBranch* Zb= tree-> Branch("Zev",&Z,"Zev/I");
  TBranch* Nb= tree-> Branch("Nev",&N,"Nev/I");
  TBranch* lb= tree-> Branch("l",&l,"l/I");


  TBranch* delEb= tree-> Branch("delE",&deltaE,"delE/F");
  TBranch* Jb= tree-> Branch("J",&J,"J/F");
  TBranch* pb= tree-> Branch("p",&p,"p/F");


  char line[500];
  char trimed_line[500];

  int i;
  //scan stdin for data about the nuclei
  while(fgets(line, 500, stdin)){
    sscanf(line, "%s", trimed_line); //removes leading whitespace
    //ignore comments '#'
    if(trimed_line[0]=='#'){
      continue;
    }  
    sscanf(line, "%d%d%f%f%d%f", &Z,&N,&J,&deltaE,&l,&p);
    //write to file and clean up
    tree->Fill();
    //printf("i: %i\n",i++);
  } 
  f.Write("",TObject::kOverwrite);
  f.Close();
}
