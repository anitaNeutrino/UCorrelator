


const char * decimated_pattern = "decimated/%d_deconv.root";
const char * wais_pattern = "wais/wais_hpol_%d_deconv.root"  ; 

int wais_start= 332;
int wais_stop = 362; 

#include "macro/cuts.C"

#include "AnitaTMVA.h" 


void doTMVA(int decimated_start = 130, int decimated_stop=439, int nworkers = 1) 
{

  // Step 1: load data

  TChain signal("wais"); 
  TChain bg("decimated"); 
  TString tmp; 

  for (int i = wais_start; i<= wais_stop; i++)
  {
    tmp.Form(wais_pattern,i); 
    signal.Add(tmp.Data()); 
  }

  for (int i = decimated_start; i<= decimated_stop; i++)
  {
    tmp.Form(decimated_pattern,i); 
    bg.Add(tmp.Data()); 
  }


  // step 1.5, enable proof if nworkers > 1

  TProof * proof = 0; 

  if (nworkers > 1) 
  {
    proof = TProof::Open("",TString::Format("workers=%d",nworkers)); 
    proof->Load("macro/proofloader.C",true); 
    proof->Load("macro/cuts.C",true); 

    signal.SetProof(); 
    bg.SetProof(); 
  }


  //Step 2: set cuts

  TCut signal_cut = isWais && isReal; 
  TCut bg_cut = thermal_sample && brightestPeak && notTooFiltered; 

  AnitaTMVA::MVAVarSet varset; 
  //setup variables
  varset.add(AnitaTMVA::MVAVar("peak.value[][]","mapPeak")); 
  varset.add(AnitaTMVA::MVAVar("peak.sigma_theta[][]","mapSigmaTheta")); 
  varset.add(AnitaTMVA::MVAVar("peak.sigma_phi[][]","mapSigmaPhi")); 
  varset.add(AnitaTMVA::MVAVar("peak.rho[][]","mapRho")); 
  varset.add(AnitaTMVA::MVAVar("coherent.peakHilbert[][]","coherentHilbertPeak")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.peakHilbert[][]","deconvHilbertPeak")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.peakHilbert[][]/totalPower[][]","deconvHilbertPeakPowerRatio")); 
  varset.add(AnitaTMVA::MVAVar("coherent_filtered.peakHilbert[][]","coherentFilteredHilbertPeak")); 
  varset.add(AnitaTMVA::MVAVar("coherent_filtered.peakHilbert[][]/coherent.peakHilbert[][]","coherentFilteredHilbertPeakRatio")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved_filtered.peakHilbert[][]","deconvFilteredHilbertPeak")); 
  varset.add(AnitaTMVA::MVAVar("coherent.peakTime[][]","coherentPeakTime")); 
  varset.add(AnitaTMVA::MVAVar("coherent.width_50_50[][]","coherentWidth5050")); 
  varset.add(AnitaTMVA::MVAVar("coherent.width_10_10[][]","coherentWidth1010")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.width_50_50[][]","deconvolvedWidth5050")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.width_10_10[][]","deconvolvedWidth1010")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.peakTime[][]","deconvolvedPeakTime")); 

  //add spectactors 
  varset.add(AnitaTMVA::MVAVar("run","run",'I',true)); 
  varset.add(AnitaTMVA::MVAVar("eventNumber","eventNumber",'I',true)); 




  //alright, this is dumb. 
  FILE * f = fopen("thermalCutTrees.root","r"); 

  if (!f) 
  {
    TFile newOut("thermalCutTrees.root","CREATE"); 
    TTree * sigtree= AnitaTMVA::makeTMVATree(&signal, &newOut,  "signal_in", varset, signal_cut); 
    TTree * bgtree= AnitaTMVA::makeTMVATree(&bg, &newOut,  "bg_in", varset, bg_cut); 
    newOut.Write(); 
  }
  else
  {
    fclose(f); 
  }



  
  TFile  out("thermalCutTrees.root"); 
  TTree* sigtree = (TTree*) out.Get("signal_in"); 
  TTree* bgtree = (TTree*) out.Get("bg_in"); 
  

  TFile tmvaOut("thermalCuts.root","RECREATE"); 
  // Step 3: setup TMVA  
  TMVA::Factory *factory = new TMVA::Factory("thermal_cuts", &tmvaOut,"V"); 


  TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 


  /// this is where the input selection cuts  are set now until I can figure out the entrylist
  
  dl->AddSignalTree(sigtree); 
  dl->AddBackgroundTree(bgtree); 

  varset.setUpData(dl); 


  //setup methods 

  factory->BookMethod(dl, TMVA::Types::kFisher, "Fisher"); 
//  factory->BookMethod(dl, TMVA::Types::kKNN, "kNN"); 
  factory->BookMethod(dl, TMVA::Types::kBDT, "BDT"); 
  factory->BookMethod(dl, TMVA::Types::kMLP, "MLP"); 


  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 

  tmvaOut.Write(); 






}


