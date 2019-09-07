

const char * data_pattern = "a3all/%d_%s.root";
const char * data_tree = "anita3"; 
const char * wais_pattern = "wais/%d_%s.root"  ; 
const char * simulated_pattern = "%s/*%s.root"; 

int wais_start= 332;
int wais_stop = 362; 


#ifndef CUTS_LOADED
#include "macro/cuts.C"
#endif
#include "AnitaTMVA.h" 


void makeTrees(int data_start = 160, int data_stop=439, const char * mc_dir = "simulated_kotera_max", const char * filter = "sinsub_10_3_ad_2", int nworkers = 8) 
{

  // Step 1: load data

  TChain signal(mc_dir ? "simulation" :"wais"); 
  std::vector<TChain*> bg(nworkers); 
  if (nworkers > 1) 
    ROOT::EnableThreadSafety(); 


  for  (int i = 0; i < nworkers; i++) 
  {
    bg[i] = new TChain(data_tree); 
  }


  TString tmp; 

  if (mc_dir) 
  {
    tmp.Form(simulated_pattern,mc_dir, filter); 
    signal.Add(tmp.Data()); 
  }
  else
  {
    for (int i = wais_start; i<= wais_stop; i++)
    {
      tmp.Form(wais_pattern,i,filter); 
      signal.Add(tmp.Data()); 
    }
  }


  for (int i = data_start; i<= data_stop; i++)
  {
    tmp.Form(data_pattern,i,filter); 
  //  printf("added %d to %d (%p)\n", i, i%nworkers, bg[i % nworkers]); 
    bg[i % nworkers]->Add(tmp.Data()); 
  }



  //Step 2: set cuts

  TCut signal_cut= mc_dir ?  isMC : isWais && isReal; //revisit this 
  TCut bg_cut = thermal_sample && brightestPeak; 

  AnitaTMVA::MVAVarSet varset("tree_vars.tmva"); 


  TString treefilename; 
  treefilename.Form("thermalCutTrees_%s_%s.root",filter,mc_dir ? mc_dir: "wais"); 

  TFile newOut(treefilename.Data(),"RECREATE"); 

  TTree * sigtree= AnitaTMVA::makeTMVATree(&signal, &newOut,  "sig_in", varset, signal_cut); 

  sigtree->Write(); 
  
  TTree * bgtree= AnitaTMVA::makeTMVATree(nworkers, (TTree**) &bg[0], &newOut, "bg_in", varset, bg_cut); 

  bgtree->Write(); 

  newOut.Write(); 
}

void doTMVA(int data_start = 160, int data_stop=439, const char * mc_dir = "simulated_kotera_max",  const char * filter = "sinsub_10_3_ad_2", int nworkers = 8) 
{

  TString treefilename; 
  treefilename.Form("thermalCutTrees_%s_%s.root",filter,mc_dir ? mc_dir : "wais"); 

  //alright, this is dumb. 
  FILE * f = fopen(treefilename.Data(),"r"); 

  if (!f) 
  {
    makeTrees(data_start, data_stop, mc_dir,filter, nworkers); 
  }
  else
  {
    fclose(f); 
  }



  TFile  out(treefilename.Data()); 
  TTree* sigtree = (TTree*) out.Get("sig_in"); 
  TTree* bgtree = (TTree*) out.Get("bg_in"); 

  TString tmvaOutName; 
  tmvaOutName.Form("thermalCuts_%s_%s.root",filter,mc_dir ? mc_dir : "wais"); 

  TFile tmvaOut(tmvaOutName.Data(),"RECREATE"); 

  TMVA::Factory *factory = new TMVA::Factory("thermal_cuts", &tmvaOut,"V"); 

  TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 

  if (mc_dir) 
  {
    dl->SetSignalWeightExpression("weight"); 
  }

 

  /* These are the variables to be used. They must have been generated already */ 
  dl->AddVariable("mapPeak"); 
  dl->AddVariable("coherentFilteredHilbertPeakRatio"); 
  dl->AddVariable("filteredFraction"); 
  dl->AddVariable("coherentHilbertPeak"); 
  dl->AddVariable("deconvHilbertPeak"); 
  dl->AddVariable("dPhiSun"); 
  dl->AddVariable("dPhiNorth"); 
  dl->AddVariable("deconvImpulsivity"); 
  dl->AddVariable("deconvLinearPolFraction"); 
  dl->AddVariable("deconvLinearPolAngle"); 
  dl->AddVariable("deconvolvedWidth1010"); 
  dl->AddVariable("deconvolvedWidth5050"); 
  dl->AddSpectator("run"); 
  dl->AddSpectator("weight"); 
  dl->AddSpectator("eventNumber"); 
  dl->AddSignalTree(sigtree); 
  dl->AddBackgroundTree(bgtree); 


  //setup methods 
  factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher","CreateMVAPdfs=true"); 
//  factory->BookMethod(dl, TMVA::Types::kKNN, "kNN"); 
  factory->BookMethod(dl, TMVA::Types::kBDT, "BDT","CreateMVAPdfs=true"); 
//  factory->BookMethod(dl, TMVA::Types::kMLP, "MLP"); 
//  factory->BookMethod(dl, TMVA::Types::kSVM, "SVM","CreateMVAPdfs=true"); 


  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 

  tmvaOut.Write(); 


}


