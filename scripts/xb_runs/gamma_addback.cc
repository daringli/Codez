#include <cstdlib>     /* exit, EXIT_FAILURE */
#include <cstdio>
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
/*
class Neighbour_list
{
public: std::vector<std::vector<short int> > neighbours;
  Neighbour_list()
  {
    char* data_filename="cb_geometry.dat";
    FILE* f=fopen(data_filename,"r");
    char line[64];
    short int i=0;
    while(fgets(line, 64, f)){
      std::vector<short int> list(7);
      list[0]=i;
      sscanf(line, "%*i %*f %*f %*c %hi %hi %hi %hi %hi %hi", &list[1],&list[2],&list[3],&list[4],&list[5],&list[6]);
      if(list[6]==163){
	list.erase(list.begin()+6);
      }
      neighbours.push_back(list);
      i++;
      //printf("i: %hi\n",i);
    }
  }
};
*/

int main(int argc, char* argv[]) 
{
  char* in_filename="rootfile.root";
  char* out_filename="gamma_addback.tsv";

  TFile f(in_filename,"update");
  float Eex;
  int Z, N;
  TTree *h102 = (TTree*)f.Get("h102");

  //read E* from a file
  char* quasi_out="qout.tsv";
  FILE* fid= fopen(quasi_out,"r");
  char line[500];

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

  int i=0;
  FILE* fout=fopen(out_filename,"w");
  fprintf(fout,"Event\tE*\tE_ab\n");

  while(fgets(line, 500, fid)){
    if(line[0]!='*'){
      continue;
    }
    //next line contains the goodies=first nuclei in each event.
    fgets(line, 500, fid);
    sscanf(line, "%i%i%f%*[^\n]",&Z,&N,&Eex); //get third entry into &E  
    //fills our addresses with value!
    int energy_bytes=Xbeb->GetEvent(i); 
    Xbtb->GetEvent(i); 
    Xbib->GetEvent(i); 
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
      puts("in while!");
      //find indices of deposits in neighbours. 
      std::vector<short int> neighbours=neighbour_list(Xbi_event[0]);
      //start at end so that we can remove without messing up future indices.
      double addback_E=0;
      for(int k=deposits_left-1;k>=0;k--){
	//if deposit k is in a neighbour and the time difference is <30ns
	//trivially true for k=0
	if((find(neighbours.begin(),neighbours.end(),Xbi_event[k])!=neighbours.end()) && (abs(Xbt_event[k]-Xbt_event[0])<30)){
	  //add the energy and remove this deposit from the vector of deposits
	  addback_E=addback_E+Xbe_event[k]; 
	  Xbe_event.erase(Xbe_event.begin()+k);
	  Xbt_event.erase(Xbt_event.begin()+k);
	  Xbi_event.erase(Xbi_event.begin()+k);
	}
      }

      addback_energies.push_back(addback_E);
      fprintf(fout,"%i\t%f\t%f\n",i+1,Eex,addback_E);
    }
    i++;
  }
  fclose(fid);
  fclose(fout);
  
  //Overwrite the old tree! I got <TBasket::Streamer>: The value of fKeylen is incorrect when I tried to create new trees. D:
  f.Write("",TObject::kOverwrite);     // save only the new version of the tree
  f.Close();
}


std::vector<short int> neighbour_list(int index){
  if(index>162 && index<0){
    //not valid indices, and thus no neighbours.
    return std::vector<short int>(0);
  }
  char data_filename[]="cb_geometry.dat";
  FILE* f=fopen(data_filename,"r");
  char line[66];
  for(int i=1;i<=index;i++){
    fgets(line, 65, f);
  }
  short int index_list[7];
  index_list[0]=(short int)index; //every crystal is its own neighbour!!!
  sscanf(line, "%*i %*f %*f %*c %hi %hi %hi %hi %hi %hi", &index_list[1],&index_list[2],&index_list[3],&index_list[4],&index_list[5],&index_list[6]); 
  int length;
  if(index_list[6]==163){
    length=6;
  }
  else{
    length=7;
  }
  std::vector<short int> ret(index_list,index_list+length);
  fclose(f);
  return ret;
}

