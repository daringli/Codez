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

#include "elements.hh"


int main(int argc, char *argv[])
{
  std::vector<std::string> Event_ranges;
  bool use_all=1; //by default, we use all the events on the list
  std::vector<int> events;

  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("events,n",  po::value<std::vector<std::string> >(&Event_ranges)->multitoken(),"Event numbers to include. Can be given as a single value or range, and given several times.") 
    ;

 
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm); 

  //if help, print description.
  if ( vm.count("help")  ) { 
    std::cout << desc << std::endl; 
    return 0; 
  } 

 if(vm.count("events")){
    use_all=0; //now we do not want to use all events, but rather the users list
    for(int i=0; i<Event_ranges.size();i++){
      int sep=Event_ranges[i].find("-");
      if(sep==std::string::npos){
	//add the single event to the list of events we want
	events.push_back(atoi(Event_ranges[i].c_str()));
      }
      else{
	//add the entire range of events to our list of events
	int start=atoi(Event_ranges[i].substr(0,sep).c_str());
	int stop=atoi(Event_ranges[i].substr(sep+1,std::string::npos).c_str());
	for(int j=start;j<=stop;j++){
	  events.push_back(j);
	}
      }
    }
  }

  //initilize event counter
  int event_counter=0;

  //get nuclei to run from stdin
  int Z,N;
  double E,J,px,py,pz;

  char line[500];
  char trimed_line[500];
  bool use_event=1;
  //scan stdin for data about the nuclei
  while(fgets(line, 500, stdin)){
    sscanf(line, "%s", trimed_line); //removes leading whitespace
    //ignore comments '#' and go to next desired event at '*'
    if(trimed_line[0]=='#'){
      continue;
    }  
    else if(trimed_line[0]=='*'){
      event_counter++; //currently ignores the actual number that is printed
      if(!use_all){
	//we must check if we want this event
	use_event=!(find(events.begin(), events.end(), event_counter)==events.end());
	  if(use_event==0){
	    //do not proceed and clean up old event if we do not want new.
	    continue;
	  }
      }
      if(event_counter>1){
	printf("}\n"); //to signify the end of the old event
      }
      printf("*** EVENT %i ***\n",event_counter);
      printf("origin: {\n");
      continue;
    }  
    //if user specified which events they want, we ignore ones they don't
    if(use_event==0){
      continue;
    }

    sscanf(line, "%d%d%lf%lf%lf%lf%lf", &Z,&N,&E,&J,&px,&py,&pz);
    Nucleus n;
    n.set_Z(Z);
    n.set_N(N);

    printf("  %s:",element_name(n));	
    printf("     pxyz=(%e, %e, %e)\n",px,py,pz);
  }
  printf("}\n"); //close final event
} 
