#include "AnitaTMVA.h" 
void runTMVA(const char * outname = "tmva/a3all_fisher_tmva_all.root", const char * signal_file = "thermalTrees/simulated_*_sinsub_10_3_ad_2.root", const char * signal_tree_name = "simulation", const char * bg_file = "thermalTrees/a3all_*_sinsub_10_3_ad_2.root", const char *  bg_tree_name = "anita3") 
{
  TChain sig(signal_tree_name); 
  sig.Add(signal_file);
  sig.SetImplicitMT(true); 

  TChain bg(bg_tree_name);
  bg.Add(bg_file); 
  bg.SetImplicitMT(true); 

  TFile tmvaOut(outname,"RECREATE"); 

  TMVA::Factory *factory = new TMVA::Factory("thermal_cuts", &tmvaOut,"V"); 

  TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 

  if (!strcmp(signal_tree_name,"simulation"))
  {
    dl->SetSignalWeightExpression("weight"); 
  }
 
  AnitaTMVA::MVAVarSet vars("tree_vars.tmva"); 
  vars.setUpData(dl); 


  dl->AddSignalTree(&sig); 
  dl->AddBackgroundTree(&bg); 
  TCut cut("!isnan(blastFraction) && !isnan(absHwAngle) && !isnan(deconvLinearPolFraction) && !isnan(absDeconvLinearPolAngle) "); 
  dl->AddCut(cut); 
  dl->AddCut(TCut("isMostImpulsive && pointsToMC && (MaxPeak < 1000) && !payloadBlast && chisq < 0.01 && ((iteration < 5 && HPolTrigger) || (iteration >=5 && VPolTrigger))"),"Signal"); 
  dl->AddCut(TCut("theta < 0 && !payloadBlast && isMostImpulsive  && chisq < 0.01 && MaxPeak < 1000 &&  ((iteration < 5 && HPolTrigger)  || (iteration >=5 && VPolTrigger))") ,"Background"); 

  //setup methods 
  factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher",""); 
//  factory->BookMethod(dl, TMVA::Types::kFisher, "TransformedFisher","VarTransform=G,D"); 
//  factory->BookMethod(dl, TMVA::Types::kKNN, "kNN"); 
//  factory->BookMethod(dl, TMVA::Types::kCuts, "Cuts"); 
//  factory->BookMethod(dl, TMVA::Types::kCuts, "Cuts","CreateMVAPdfs=true"); 
//  factory->BookMethod(dl, TMVA::Types::kBDT, "BDT"); 
//  factory->BookMethod(dl, TMVA::Types::kMLP, "MLP"); 
//  factory->BookMethod(dl, TMVA::Types::kSVM, "SVM"); 


  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 
  tmvaOut.Write(); 
}


