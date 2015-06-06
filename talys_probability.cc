#include <cstdio>
#include <string>
#include <string.h> //string
#include <fstream>
#include <iostream>
#include <sstream>

#include "decay.hh"
#include "decay_name.hh"
#include "talys_probability.hh"


std::vector<Particle_probability> talys_width(std::string talys_outfile_name)
{
  //read a Talys outfile and parse the first table in it
  FILE* talys_outfile = fopen(talys_outfile_name.c_str(),"r");

  //initialize widths
  std::vector<Particle_probability> ret(7);
  for(int i=0;i<7;i++){
    ret[i].prob=0;
  }

  int i=1;
  char line[500];

  //loop through every line in file
  while(fgets(line, 500, talys_outfile)){
    //first 5 lines are header lines
    if(i<6){
      //on line 4 is the channel labels
      if(i==4){
	//read the channel labels with stringstream magic!
	std::string readto;
	std::string linestring=std::string(line);
	std::stringstream linestream(linestring);
	int j=1;
	while(linestream>>readto)
	  {
	    //printf("%s\n",readto.c_str());
	    //skip first header column, since it doesn't name channel
	    if(j>1){
	      if(readto=="proton"){
		readto="p";
	      }
	      else if(readto=="neutron"){
		readto="n";
	      }
	      else if(readto=="deuteron"){
		readto="d";
	      }
	      else if(readto=="triton"){
		readto="t";
	      }
	      else if(readto=="helium-3"){
		readto="He3";
	      }
	      ret[j-2].name=readto;
	    }
	    j++;
	  }
      }
      i++;
      continue;
    }
    //after 5th line, we read each line with stringstream magic!
    int j=1;
    double fnumber; //placeholder float to read partial widths into
    std::string linestring=std::string(line);
    std::stringstream linestream(linestring);
    while(linestream>>fnumber){
      //skip first column, since energy rather than width
      if(j>1){
	ret[j-2].prob=ret[j-2].prob+fnumber;
	//printf("fnumber: %f\n",fnumber);
      }
      j++;
    }
    if(j==1){
      goto endloop; //the table ended, and thus we stop reading!
    }
    i++;
  }
 endloop:
  //normalize decay probabilities
  double probtot=0;
  for(int j=0; j<7;j++){
    probtot=probtot + ret[j].prob;
  }
  for(int j=0; j<7;j++){
    ret[j].prob=ret[j].prob/probtot;
    //printf("branch %s: %f\n",ret[j].name.c_str(),ret[j].prob);
  }
  fclose(talys_outfile);
  return ret;
}

