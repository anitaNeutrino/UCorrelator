
const char * data_tree = "anita4"; 
const char * data_pattern = "/Volumes/SDCard/data/a4all/%d_%s.root";
const char * wais_pattern = "/Volumes/SDCard/data/wais/%d_%s.root"  ; 
const char * simulated_pattern = "/Volumes/SDCard/data/simulated/%d_%s.root"; 

int wais_start= 120;
int wais_stop = 155; 

int simulation_start= 1;
int simulation_stop = 200; 


#ifndef CUTS_LOADED
#include "macro/cuts.C"
#endif
#include "AnitaTMVA.h" 


void makeTrees(int data_start = 50, int data_stop=367, const char * mc_dir = "simulated", const char * filter = "max_30002_sinsub_10_3_ad_2", int nworkers = 8) 
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
    for (int i = simulation_start; i<= simulation_stop; i++)
    {
      tmp.Form(simulated_pattern,i, "max_1000_sinsub_10_3_ad_2"); 
      signal.Add(tmp.Data()); 
    }
    std::cout<<"simulation entries: "<< signal.GetEntries() << std::endl;
    
  }
  else
  {
    for (int i = wais_start; i<= wais_stop; i++)
    {
      tmp.Form(wais_pattern,i,filter); 
      signal.Add(tmp.Data()); 
    }
    std::cout<<"wais entries: "<< signal.GetEntries() << std::endl;
  }


  for (int i = data_start; i<= data_stop; i++)
  {
    tmp.Form(data_pattern,i,filter); 
  //  printf("added %d to %d (%p)\n", i, i%nworkers, bg[i % nworkers]); 
    bg[i % nworkers]->Add(tmp.Data()); 
  }
  std::cout<<"background entries: "<< bg[0]->GetEntries() << std::endl;



  //Step 2: set cuts

  TCut signal_cut= mc_dir ?  mc_sample : wais_sample; //revisit this 
  TCut bg_cut = thermal_sample && fromSky; 

  AnitaTMVA::MVAVarSet varset("tree_vars.tmva"); 


  TString treefilename; 
  treefilename.Form("thermalCut/thermalCutTrees_%s_%s.root",filter,mc_dir ? mc_dir: "wais"); 

  TFile newOut(treefilename.Data(),"RECREATE"); 

  TTree * sigtree= AnitaTMVA::makeTMVATree(&signal, &newOut,  "sig_in", varset, signal_cut); 

  sigtree->Write(); 
  
  TTree * bgtree= AnitaTMVA::makeTMVATree(nworkers, (TTree**) &bg[0], &newOut, "bg_in", varset, bg_cut); 

  bgtree->Write(); 

  newOut.Write(); 
}

void doTMVA(int data_start = 50, int data_stop=367, const char * mc_dir = "simulated",  const char * filter = "max_30001_sinsub_10_3_ad_2", int nworkers = 8) 
{

  TString treefilename; 
  treefilename.Form("thermalCut/thermalCutTrees_%s_%s.root",filter,mc_dir ? mc_dir : "wais"); 

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
  tmvaOutName.Form("thermalCut/thermalCuts_%s_%s.root",filter,mc_dir ? mc_dir : "wais"); 

  TFile tmvaOut(tmvaOutName.Data(),"RECREATE"); 

  TMVA::Factory *factory = new TMVA::Factory("thermal_cuts", &tmvaOut,"V"); 
  TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 
  if (mc_dir) 
  {
    dl->SetSignalWeightExpression("weight"); 
  }
   /* These are the variables to be used. They must have been generated already */ 
  dl->AddVariable("mapPeak"); 
  dl->AddVariable("mapSNR"); 
  dl->AddVariable("deltaHilbertPeak"); 
  dl->AddVariable("deconvHilbertPeak"); 
  dl->AddVariable("deconvImpulsivity"); 
  dl->AddVariable("deconvLinearPolFraction"); 
  // dl->AddVariable("secondPeakRatio");
  dl->AddSpectator("theta"); 
  dl->AddSpectator("run"); 
  dl->AddSpectator("weight"); 
  dl->AddSpectator("eventNumber"); 
  dl->AddSpectator("evnum1"); 
  dl->AddSpectator("evnum2"); 
  dl->AddSpectator("pol"); 
  dl->AddSpectator("peak"); 
  // dl->AddSpectator("countChan"); 
  dl->AddSignalTree(sigtree,1.0); 
  dl->AddBackgroundTree(bgtree,1.0);
  //setup methods 
  // factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher","CreateMVAPdfs=true");
  // factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher0","VarTransform=N");
  // factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher1","VarTransform=D+N");
  // factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher2","VarTransform=D+G+N");
  // factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher3","VarTransform=D+G+D+G+N");

  factory->BookMethod(dl, TMVA::Types::kLD, "LD0","CreateMVAPdfs=true");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD0","VarTransform=N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD1","VarTransform=D+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD2","VarTransform=D+G+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD4","VarTransform=G+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD5","VarTransform=G+D+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD6","VarTransform=G+D+G+D+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD7","VarTransform=P+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD8","VarTransform=P+D+N");
  // factory->BookMethod(dl, TMVA::Types::kLD, "LD9","VarTransform=D_Signal+G_Signal+N_Signal");

// factory->BookMethod(dl, TMVA::Types::kMLP, "MLPBFGS","H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS" ); 
 // factory->BookMethod(dl, TMVA::Types::kKNN, "kNN"); 
// factory->BookMethod(dl, TMVA::Types::kBDT, "BDT","CreateMVAPdfs=true"); 
//  factory->BookMethod(dl, TMVA::Types::kMLP, "MLP"); 
 // factory->BookMethod(dl, TMVA::Types::kSVM, "SVM","CreateMVAPdfs=true"); 


  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 



  tmvaOut.Write(); 


}


