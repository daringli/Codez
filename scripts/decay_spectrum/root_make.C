//cint script, run by .L root_make.C
{
  TTree* thetree=new TTree();
  thetree->ReadFile("../../test.tsv","Z/I:N/I:l/I:Jf/F:Ef/F");
  TCanvas* theCanvas=new TCanvas("h1","h1",600,600);
  thetree->Draw("Ef","Z==2 && N==2"); //final energy after alpha decay
}
