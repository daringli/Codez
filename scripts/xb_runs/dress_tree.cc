#include <cstdlib>     /* exit, EXIT_FAILURE */
#include <cstdio>
#include <cstring>
#include <vector>
#include <cstring>
#include <algorithm>

#include "TColor.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TTree.h"
#include "TSpectrum.h"
#include "THStack.h"
#include "TLegend.h"
#include "TFriendElement.h"
#include "TStyle.h"
#include "TPaveLabel.h"

std::vector<short int> neighbour_list(int index);

int main() 
{
  //makes the tree nice and shiny
  //with prefragment E* and addbacked Xbe!
  char* in_filename="rootfile.root";

  TFile f(in_filename,"update");
  float Eex;
  int Z, N;
  TTree *h102 = (TTree*)f.Get("h102;1");



  float Xbe[50]; 
  float Xbt[50]; //addresses to store energy and time
  int Xbi[50]; //address to store crystal index
  TBranch *Xbeb = h102->GetBranch("Xbe");
  TBranch *Xbtb = h102->GetBranch("Xbt");
  TBranch *Xbib = h102->GetBranch("Xbi");

  Xbeb->SetAddress(Xbe);
  Xbtb->SetAddress(Xbt);
  Xbib->SetAddress(Xbi);

  std::vector<float> addback_energies(0);
  addback_energies.reserve(50); //so that our address used below does not get invalidated by push_back.
  int num_addback; //the length of our addback energies above. To be assigned later.
  int fg; //fired
  //This is where we want to add new values.
  TBranch* numaddsb=h102-> Branch("num_addback",&num_addback,"num_addback/I"); //DEALBREAKER?
  TBranch* addeb= h102-> Branch("e_addback",&addback_energies[0],"addback_e[num_addback]/F");
  TBranch* eexb= h102-> Branch("Eex",&Eex,"Eex/F");
  TBranch* fgb= h102-> Branch("fg",&fg,"fg/I");


  int i=0;
  int n_events=h102->GetEntries();
  printf("Entries before: %i\n",n_events);

  //first a hackish "function" to get the number of gammas fired in each event
  std::vector<int> n_gammas(0);
 
  int event_nr=-1;
  FILE* ffgun=fopen("test.gun","r");
  char line[500];

  while(fgets(line, 500, ffgun)){
    if(line[0]=='*'){
      event_nr++;
      n_gammas.push_back(0);
    }
    if(strlen(line)>=3){
      if(line[2]=='g'){
	n_gammas[event_nr]=n_gammas[event_nr]+1;
      }
    }
  }
  fclose(ffgun);
  

  FILE* out_fg=fopen("out_fg.tsv","w");
  for(int dumi=0;dumi<event_nr;dumi++){
    fprintf(out_fg,"%i\n",n_gammas[dumi]);
  }
  fclose(out_fg);
 
  //read E* from a file
  char* quasi_out="qout.tsv";
  FILE* fid= fopen(quasi_out,"r");
  //char line[500];
  while(fgets(line, 500, fid)){
    //printf("total loops: %i\n",iiiii++);
    //printf("event loops: %i\n",i);


    if(line[0]!='*'){
      continue; //excitation energy of this event is on another line
    }
    //next line contains the goodies=first nuclei in each event.
    fgets(line, 500, fid);
    sscanf(line, "%i%i%f",&Z,&N,&Eex); //get third entry into &E  
    int energy_bytes=Xbeb->GetEvent(i); //fills our addresses with value!
    Xbtb->GetEvent(i); //fills our addresses with value!
    Xbib->GetEvent(i); //fills our addresses with value!
    int num_deposits=energy_bytes/sizeof(1.0f); //assumed to be same size
    //read the energy deposits of this event into a std::vector!
    std::vector<float> Xbe_event(Xbe, Xbe + num_deposits);
    std::vector<float> Xbt_event(Xbt, Xbt + num_deposits);
    std::vector<int> Xbi_event(Xbi, Xbi + num_deposits);
    //this list of energy deposits is already sorted so that highest e is first.
    
    //loop through deposits from highest to lowest energy, and add
    //energy in nearest neighbours if the time difference is low.
    //note: crystal its own neighbour, and hence deposit under
    //consideration will always be added.
    //included deposits are then removed from the vector, which both shrinks it
    //and changes the [0] element.
    int deposits_left;
    while((deposits_left=Xbi_event.size())>0){
      //always use zero, since we update list by removing previously used deposits.
      std::vector<short int> neighbours=neighbour_list(Xbi_event[0]); 
      //find indices of deposits in neighbours. 
      //start at end so that we can remove without messing up future indices.
      double addback_E=0;
      for(int k=deposits_left-1;k>=0;k--){
	//if deposit k is in a neighbour and the time difference is <30ns
	//trivially true for k=0
	if((find(neighbours.begin(),neighbours.end(),Xbi_event[k])!=neighbours.end()) && (abs(Xbt_event[k]-Xbt_event[0])<30)){
	  //if the energy in k is under 10MeV...
	  if(Xbe_event[k]>10 && Xbe_event[k]<0.1){
	    Xbe_event.erase(Xbe_event.begin()+k);
	    Xbt_event.erase(Xbt_event.begin()+k);
	    Xbi_event.erase(Xbi_event.begin()+k);
	    continue;
	  }
	  //...add the energy and remove this deposit from the vector of deposits!
	  //(if not, it got removed)

	  addback_E=addback_E+Xbe_event[k];
	  Xbe_event.erase(Xbe_event.begin()+k);
	  Xbt_event.erase(Xbt_event.begin()+k);
	  Xbi_event.erase(Xbi_event.begin()+k);
	}
      }
      addback_energies.push_back(addback_E);
    }
    num_addback=addback_energies.size();
    fg=n_gammas[i];
    //add the addback energies to the tree as this event.
    //fill new branches of the tree
    numaddsb->Fill(); 
    addeb->Fill(); 
    eexb->Fill(); 
    fgb->Fill(); 

    addback_energies.resize(0); //clean up for next event, but hopefully keep address to first element valid (keeping capacity is thus important, too).
    //addback_energies.reserve(50);
    //purge Eex too?
  
  i++;
  }
fclose(fid);
printf("Entries after: %i\n",h102->GetEntries());

//Overwrite the old tree! I got <TBasket::Streamer>: The value of fKeylen is incorrect when I tried to create new trees. D:
 f.Write("", TObject::kOverwrite);
f.Close();

}


std::vector<short int> neighbour_list(int index){
  if(index>162 && index<0){
    //not valid indices, and thus no neighbours.
    return std::vector<short int>(0);
  }
  char data_filename[]="cb_geometry.dat";
  FILE* f=fopen(data_filename,"r");
  char line[64];
  for(int i=1;i<=index;i++){
    fgets(line, 64, f);
  }
  int index_list[7];
  index_list[0]=(short int)index; //every crystal is its own neighbour!!!
  sscanf(line, "%*i %*f %*f %*c %hi %hi %hi %hi %hi %hi", index_list+1,index_list+2,index_list+3,index_list+4,index_list+5,index_list+6); 
  int length;
  if(index_list[5]==163){
    length=6;
  }
  else{
    length=7;
  }
  fclose(f);
  std::vector<short int> ret(index_list,index_list+length);
  return ret;
}



/*

  int main(int argc, char* argv[]){
  
  //-----read and merge root files----
  char* in_filename="rootfile.root";
  TFile *in_tfile = TFile::Open(in_filename);
  //get TTree from a rootfile by giving its label.
  tree=(TTree*)in_tfile->Get("h102");
  
  //create a friend tree to fill with Eex.
  TTree* f_tree = new TTree("T","Friend");


 

  //-----Add a new leaf to friend with entries for each event!-----
  double E;

  TBranch *beex = f_tree->Branch("Eex",&E,"Eex/F");
  //merged_tree->SetBranchAddress("Eex",&E);


  char* quasi_out="qout.tsv";
  FILE* fid= fopen(quasi_out,"r");
  char line[500];
  char trimed_line[500];

  int Z, N;
  int i=1;
  //scan lines for the energy of the first nuclei in each event.
  while(fgets(line, 500, fid)){
  sscanf(line, "%s", trimed_line); //removes leading whitespace
  //ignore comments '#' and go to next event at '*'
  if(trimed_line[0]=='*'){
  //next line contains the goodies=first nuclei in each event.
  fgets(line, 500, fid);
  sscanf(line, "%i%i%f",&Z,&N,&E); //get third entry into &E  
  f_tree->Fill(); //fill tree with content of *(&E)
  printf("i:%i\n",i);
  i++;
  }  
  }
  //----Write output to a new root file------------
  //name rootfile <in_filename>+"2"
  char out_filename[50];
  strcpy(out_filename,in_filename);
  strcat(out_filename,"2");
  TFile *f = new TFile(out_filename,"RECREATE");
  


  //merged_tree->Print();
  tree->Write();
  f->Write();
  delete f;
  delete tree;
  }
*/
