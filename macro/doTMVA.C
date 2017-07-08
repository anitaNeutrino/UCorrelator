


const char * decimated_pattern = "decimated/%d_%s.root";
const char * wais_pattern = "wais/%d_%s.root"  ; 
const char * simulated_pattern = "simulated/%d_%s.root"; 

int wais_start= 332;
int wais_stop = 362; 


#include "macro/cuts.C"
#include "AnitaTMVA.h" 


void makeTrees(int decimated_start = 130, int decimated_stop=439, int mc_run = 223, const char * filter = "sinsub_10_3_ad_2", int nworkers = 1) 
{

  // Step 1: load data

  TChain signal(mc_run ? "simulation" :"wais"); 
  TChain bg("decimated"); 
  TString tmp; 

  if (mc_run) 
  {
    tmp.Form(simulated_pattern,mc_run, filter); 
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


  for (int i = decimated_start; i<= decimated_stop; i++)
  {
    tmp.Form(decimated_pattern,i,filter); 
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

  TCut signal_cut= mc_run ?  isMC && brightestPeak : isWais && isReal; //revisit this 
  TCut bg_cut = thermal_sample && brightestPeak; 

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
  varset.add(AnitaTMVA::MVAVar("deconvolved.width_50_50[][]","deconvolvedWidth5050")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.width_10_10[][]","deconvolvedWidth1010")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.peakTime[][]","deconvolvedPeakTime")); 
  varset.add(AnitaTMVA::MVAVar("deconvolved.impulsivityMeasure[][]","deconvImpulsivity")); 
  varset.add(AnitaTMVA::MVAVar("sqrt(TMath::Power(deconvolved[][].Q,2) + TMath::Power(deconvolved[][].U,2))/deconvolved[][].I","deconvLinearPolFraction")); 
  varset.add(AnitaTMVA::MVAVar("TMath::ATan2(deconvolved[][].U,deconvolved[][].Q)/90/TMath::Pi()","deconvLinearPolAngle")); 
  varset.add(AnitaTMVA::MVAVar("abs(FFTtools::wrap(peak.phi-sun.phi,360,0))","dPhiSun")); 
  varset.add(AnitaTMVA::MVAVar("abs(FFTtools::wrap(peak.phi-heading,360,0))","dPhiNorth")); 
  varset.add(AnitaTMVA::MVAVar("flags.meanPowerFiltered[0] / flags.meanPower[0]","filteredFraction")); 

  //add spectactors . Actually here it doesn't matter if they're spectators or not I don't think 
  varset.add(AnitaTMVA::MVAVar("run","run",'I',true)); 
  varset.add(AnitaTMVA::MVAVar("eventNumber","eventNumber",'I',true)); 

  TString treefilename; 
  treefilename.Form("thermalCutTrees_%s_mc%d.root",filter,mc_run); 


  TFile newOut(treefilename.Data(),"RECREATE"); 
  TTree * sigtree= AnitaTMVA::makeTMVATree(&signal, &newOut,  "sig_in", varset, signal_cut); 
  TTree * bgtree=  AnitaTMVA::makeTMVATree(&bg, &newOut,  "bg_in", varset, bg_cut); 
  newOut.Write(); 
}

void doTMVA(int decimated_start = 390, int decimated_stop=420, int mc_run = 19,  const char * filter = "sinsub_10_3_ad_2", int nworkers = 1) 
{


  TString treefilename; 
  treefilename.Form("thermalCutTrees_%s_mc%d.root",filter,mc_run); 

  //alright, this is dumb. 
  FILE * f = fopen(treefilename.Data(),"r"); 

  if (!f) 
  {
    makeTrees(decimated_start, decimated_stop, mc_run,filter, nworkers); 
  }
  else
  {
    fclose(f); 
  }



  TFile  out(treefilename.Data()); 
  TTree* sigtree = (TTree*) out.Get("sig_in"); 
  TTree* bgtree = (TTree*) out.Get("bg_in"); 

  TString tmvaOutName; 
  tmvaOutName.Form("thermalCuts_%s_mc%d.root",filter,mc_run); 

  TFile tmvaOut(tmvaOutName.Data(),"RECREATE"); 

  TMVA::Factory *factory = new TMVA::Factory("thermal_cuts", &tmvaOut,"V"); 

  TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 
  

  /* These are the variables to be used. They must have been generated already */ 
  dl->AddVariable("mapPeak"); 
  dl->AddVariable("coherentHilbertPeak"); 
  dl->AddVariable("deconvHilbertPeak"); 
  dl->AddVariable("dPhiSun"); 
  dl->AddVariable("deconvImpulsivity"); 
  dl->AddVariable("deconvLinearPolFraction"); 
  dl->AddVariable("deconvLinearPolAngle"); 
  dl->AddVariable("deconvolvedWidth1010"); 
  dl->AddSpectator("run"); 
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


