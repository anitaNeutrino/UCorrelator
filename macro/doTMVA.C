

const char * data_pattern = "a4all/%d_%s.root";
const char * data_tree = "anita4"; 
const char * wais_pattern = "wais/%d_%s.root"  ; 
const char * simulated_pattern = "%s/*%s.root"; 

int wais_start= 120;
int wais_stop = 155; 


#ifndef CUTS_LOADED
#include "macro/cuts.C"
#endif
#include "AnitaTMVA.h" 


void makeTrees(int data_start = 50, int data_stop=367, const char * mc_dir = "simulated_kotera_max", const char * filter = "max_30001_sinsub_10_3_ad_2", int nworkers = 8) 
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

void doTMVA(int data_start = 50, int data_stop=367, const char * mc_dir = "simulated_kotera_max",  const char * filter = "max_30001_sinsub_10_3_ad_2", int nworkers = 8) 
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

  if(ROOT_VERSION_CODE >= ROOT_VERSION(6,8,0)){
    // TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 
    if (mc_dir) 
    {
      dl->SetSignalWeightExpression("weight"); 
    }
    /* These are the variables to be used. They must have been generated already */ 
    // dl->AddVariable("mapPeak");
    dl->AddVariable("mapSNR"); 
    dl->AddVariable("coherentHilbertPeak"); 
    dl->AddVariable("deconvHilbertPeak"); 
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
  }else{
    if (mc_dir) 
    {
      factory->SetSignalWeightExpression("weight"); 
    }
     /* These are the variables to be used. They must have been generated already */ 
    // factory->AddVariable("mapPeak"); 
    factory->AddVariable("mapSNR"); 
    factory->AddVariable("deconvHilbertPeak"); 
    factory->AddVariable("deconvImpulsivity"); 
    factory->AddVariable("deconvLinearPolFraction"); 
    factory->AddVariable("secondPeakRatio");
    factory->AddSpectator("theta"); 
    factory->AddSpectator("run"); 
    factory->AddSpectator("weight"); 
    factory->AddSpectator("eventNumber"); 
    factory->AddSignalTree(sigtree); 
    factory->AddBackgroundTree(bgtree);
  }

 

  //setup methods 
  factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher","CreateMVAPdfs=true");
// factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS","H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS" ); 
//  factory->BookMethod(dl, TMVA::Types::kKNN, "kNN"); 
// factory->BookMethod(dl, TMVA::Types::kBDT, "BDT","CreateMVAPdfs=true"); 
//  factory->BookMethod(dl, TMVA::Types::kMLP, "MLP"); 
//  factory->BookMethod(dl, TMVA::Types::kSVM, "SVM","CreateMVAPdfs=true"); 


  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 

  tmvaOut.Write(); 


}


