#include <stdlib.h>     /* exit, EXIT_FAILURE */

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

int main(int argc, char* argv[]){

  if(argc!=2){
    printf("Error, incorrect number of arguments. Usage:\n %s filename\n",argv[0]);
    exit(1);
  }
  char* filename=argv[1];
  //char* drawstring=argv[2];
  //char* cuts=argv[3];


  TTree* thetree=new TTree();
  thetree->ReadFile(filename,"Z/I:N/I:l/I:Jf/F:Ef/F");

  //name rootfile filename.root
  char rootfile[50];
  strcpy(rootfile,filename);
  strcat(rootfile,".root");

  TFile *f = new TFile(rootfile,"UPDATE");
  //write a tree with the label filename to the file test.root, 
  //i.e. filename does not refer to the rootfile, 
  //but to the original file, and is just a label.
  thetree->Write(filename); 

  //TCanvas* theCanvas=new TCanvas("h1","h1",600,600);
  //thetree->Draw(drawstring,cuts);

  //TH1F* h1 = (TH1F*)gDirectory->Get("h1");
  //TAxis* ax=h1->GetXaxis();
  //TAxis* ay=h1->GetYaxis();
  //h1->SetTitle("");  
}
