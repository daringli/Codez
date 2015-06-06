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
  std::vector<std::string> Event_ranges;
  std::string out_filename;
  bool use_all=1; //by default, we use all the events on the list
  std::vector<int> events;

  namespace po = boost::program_options; 
  po::options_description desc("Options"); 
  desc.add_options() 
    ("help,h", "Print help messages") 
    ("filename,f",  po::value<std::string> (&out_filename),"Filename to write to. ") 
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

  if(vm.count("filename")==0){
    puts("Need a filename to write to!");
    return 1;
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



  //prepare rootfile
  TFile f(out_filename.c_str(),"recreate");


  //get nuclei to run from stdin. At most max_parts each event!
  int max_parts=50;
  int Z[max_parts],N[max_parts];
  float E,J,px[max_parts],py[max_parts],pz[max_parts]; 

  int part_counter=0;

  //prepare the tree and its branches
  TTree *tree = new TTree("eg_out","Event generator output");
  TBranch* num_partb=tree-> Branch("num_part",&part_counter,"num_part/I");
  TBranch* Zb= tree-> Branch("Z",&Z[0],"Z[num_part]/I");
  TBranch* Nb= tree-> Branch("N",&N[0],"N[num_part]/I");
  TBranch* pxb= tree-> Branch("px",&px[0],"px[num_part]/F");
  TBranch* pyb= tree-> Branch("py",&py[0],"py[num_part]/F");
  TBranch* pzb= tree-> Branch("pz",&pz[0],"pz[num_part]/F");


  char line[500];
  char trimed_line[500];
  bool use_event=1;

  //initilize counters
  int event_counter=0;
  //scan stdin for data about the nuclei
  while(fgets(line, 500, stdin)){
    sscanf(line, "%s", trimed_line); //removes leading whitespace
    //ignore comments '#' and go to next desired event at '*'
    if(trimed_line[0]=='#'){
      continue;
    }  
    else if(trimed_line[0]=='*'){
      if(event_counter>0){
	tree->Fill();
	part_counter=0;
      }

      event_counter++; //currently ignores the actual number that is printed
      if(!use_all){
	//we must check if we want this event
	use_event=!(find(events.begin(), events.end(), event_counter)==events.end());
	if(use_event==0){
	  //do not proceed and clean up old event if we do not want new.
	  continue;
	}
      }
      continue;
    }
    //if user specified which events they want, we ignore ones they don't
    if(use_event==0){
      continue;
    }

    sscanf(line, "%d%d%f%f%f%f%f", &Z[part_counter],&N[part_counter],&E,&J,&px[part_counter],&py[part_counter],&pz[part_counter]);
    //printf("Z[%i] = %i\n",part_counter,Z[part_counter]);
    //printf("px[%i] = %f\n",part_counter,px[part_counter]);

    part_counter++;
    if(part_counter>=max_parts){
      printf("\nToo many particles generated. Change max_part in out2root and try again!\n");
      return 1;
    }
  }
  //write to file and clean up
  tree->Fill();
  f.Write("",TObject::kOverwrite);
  f.Close();
} 
