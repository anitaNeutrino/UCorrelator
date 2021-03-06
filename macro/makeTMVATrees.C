
#ifndef CUTS_LOADED
#include "macro/cuts.C"
#endif

#include "AnitaTMVA.h" 



void makeTMVATrees(const char * dir, const char * tree_name, int start_run, int end_run, int nworkers = 1, int runs_per_file=40, const char * filter = "sinsub_10_3_ad_2") 
{

  int current_run = start_run; 

  while (current_run <= end_run) 
  {

    int end = current_run + runs_per_file-1; 

    std::vector<TChain *> c(nworkers); 

    if (nworkers > 1) 
    {
      ROOT::EnableThreadSafety(); 
    }

    for (int i = 0; i < nworkers; i++)
    {
      c[i] = new TChain(tree_name); 
    }


    TString tmp; 
    for (int i = current_run; i <= end; i++)
    {
      tmp.Form("%s/%d_%s.root", dir, i, filter); 
      c[i % nworkers]->Add(tmp.Data()); 
    }

    TCut cut = !strcmp(tree_name,"wais")       ? isWais && isReal : isReal ;

    AnitaTMVA::MVAVarSet varset("tree_vars.tmva"); 

    TString treefilename; 
    TString sanitized_dir(dir); 
    sanitized_dir.ReplaceAll("/","_"); 
    treefilename.Form("thermalTrees/%s_%d-%d_%s.root",sanitized_dir.Data(),current_run,end, filter); 

    TFile newOut(treefilename.Data(),"RECREATE"); 
    TTree * tree = AnitaTMVA::makeTMVATree(nworkers, (TTree**) &c[0], &newOut, tree_name, varset, cut); 

    tree->Write(); 
    newOut.Write(); 

    current_run += runs_per_file; 
  }
}
