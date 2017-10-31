
//script of using TMVA/MLP regression method to find the deltaT between each pair of antennas. 
//MLP is neutron network. Works like a black box to fit any arbitrary funciton. Here theta and phi are inputs, deltaT is output. 
//The network is really small and fits very well.
// //Author: Peng Cao
//needs a different correlaiton class "AllCorrelationSummaryAnita4.h" in peng's branch
const char * corr_pattern = "corrTrees/run%dCorrTree_0.root"; 

#ifndef CUTS_LOADED
#include "../macro/cuts.C"
#endif
#include "AnitaTMVA.h"
#include "AllCorrelationSummaryAnita4.h"
void doMLP_pair(int data_start = 135, int data_stop=152,int antenna1=0, int antenna2=0, int polariztion=0) 
{
  //Step 1: load data
  TString tmp; 
  TChain* chain = new TChain("corrTree"); 
  for (int i = data_start; i<= data_stop; i++){
    std::cout<< "load run "<< i <<std::endl;
    tmp.Form(corr_pattern,i); 
    chain->Add(tmp.Data());
  }
  TFile *file = TFile::Open("corrTrees/allWais.root","RECREATE");
  chain->CloneTree(-1,"fast");
  file->Write();

  float combo1 = 0;
  TTree *corrTree = (TTree*)file->Get("corrTree");
  AllCorrelationSummaryAnita4 * corr = 0;
  int old_pol;
  corrTree->SetBranchAddress("corr", &corr);
  corrTree->SetBranchAddress("pol", &old_pol);
  
  // corrTree->Print();
  corrTree->Write();

   TFile f("TMVA/pairCorTree.root","recreate");
   TTree pairCorTree("pairCorTree","a Tree with data from");


   int ant1,ant2, isValid, pol;
   double maxCorVals,maxCorTimes,expectedDeltaT,phiWave,thetaWave;
   int run, eventNumber,realTime;
   pairCorTree.Branch("ant1",&ant1);
   pairCorTree.Branch("ant2",&ant2);
   pairCorTree.Branch("isValid",&isValid);
   pairCorTree.Branch("pol",&pol);
   pairCorTree.Branch("maxCorTimes",&maxCorTimes);
   pairCorTree.Branch("maxCorVals",&maxCorVals);
   pairCorTree.Branch("expectedDeltaT",&expectedDeltaT);
   pairCorTree.Branch("phiWave",&phiWave);
   pairCorTree.Branch("thetaWave",&thetaWave);
   // pairCorTree.Branch("run",&run);
   pairCorTree.Branch("eventNumber",&eventNumber);
   // pairCorTree.Branch("realTime",&realTime);


   Long64_t maxEntry = corrTree->GetEntries();
  for(Long64_t entry=0; entry < maxEntry; entry++)
  {
    corrTree->GetEntry(entry);
    for (int i = 0; i < 78; i++){
      for(int polInd=0; polInd < 2; polInd++){
        ant1 = corr->firstAnt[i];
        ant2 = corr->secondAnt[i];  
        if(ant1 > ant2){
          int tmp = ant1;
          ant1 = ant2;
          ant2 = tmp;
        } 

        pol = 1 - old_pol;
        maxCorTimes = corr->maxCorTimes[i];
        maxCorVals = corr->maxCorVals[i];
        expectedDeltaT = corr->expectedDeltaT[i];
        phiWave = corr->phiWave;
        thetaWave = corr->thetaWave;
        if (std::abs(expectedDeltaT - maxCorTimes)<1 && std::abs(FFTtools::wrap((phiWave - corr->fAntPhi[i][0])*180.0/3.14, 360, 0))<45 && std::abs(FFTtools::wrap((phiWave - corr->fAntPhi[i][1])*180.0/3.14, 360, 0))<45){
          isValid = 1;
        }else{
          isValid = 0;
        }

        // run = corr->run;
        eventNumber = corr->eventNumber;
        // realTime = corr->realTime;
        pairCorTree.Fill();
      }
    }
  }
  pairCorTree.Write();


// // load data directly from disk
//   TString tmp = "TMVA/pairCorTree.root"; 
//   TChain pairCorTree("pairCorTree"); 
//   pairCorTree.Add(tmp.Data());



  //Step 2: output file
  TString tmvaOutName; 
  tmvaOutName.Form("TMVA/antPair_%d_%d.root",antenna1,antenna2); 
  TFile tmvaOut(tmvaOutName.Data(),"RECREATE"); 

  char weightFile[100];
  sprintf(weightFile, "WeightFile_%d_%d_%d", antenna1, antenna2, polariztion);
  (TMVA::gConfig().GetIONames()).fWeightFileExtension = weightFile;
  TMVA::Factory *factory = new TMVA::Factory("trainTiming", &tmvaOut,"!V"); 
  factory->AddRegressionTree(&pairCorTree); // corrTree is the Tchain contains all the readed trees.
  // factory->PrepareTrainingAndTestTree( "pol==0 && isValid==1 && ant1==0 && ant2==17", "SplitMode=random:!V" );
  char cutString[100];
  sprintf (cutString, "isValid==1 && ant1==%d && ant2==%d && pol==%d", antenna1, antenna2, polariztion);
  factory->PrepareTrainingAndTestTree( cutString, "SplitMode=random:!V" );

  // TMVA::DataLoader *dl = new TMVA::DataLoader("thermal"); 

  /* These are the variables to be used. They must have been generated already */ 
  // factory->AddVariable("ant1","","",'I');
  // factory->AddVariable("ant2","","",'I');
  factory->AddSpectator("ant1");
  factory->AddSpectator("ant2");
  factory->AddSpectator("isValid");
  factory->AddSpectator("pol");
  factory->AddTarget("maxCorTimes"); 
  factory->AddSpectator( "maxCorVals" );
  factory->AddSpectator( "expectedDeltaT" );
  factory->AddVariable("thetaWave"); 
  factory->AddVariable("phiWave"); 
  // factory->AddSpectator("firstAnt[0]");
  // factory->AddSpectator("secondAnt[0]"); 
  // factory->AddSpectator("combo1"); 
  // factory->AddSpectator("centreAntenna"); 

  // factory->AddSpectator("pol"); 
  // factory->AddSpectator("fAntPhi"); 


// Artificial Neural Network (Multilayer perceptron) - TMVA version
// factory->BookMethod( TMVA::Types::kMLP, "MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5" );
// NN with BFGS quadratic minimisation
// factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS" );
// Artificial Neural Network (Multilayer perceptron) - TMVA version
// factory->BookMethod( TMVA::Types::kMLP, "MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5" );
// NN with BFGS quadratic minimisation
factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS","H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N+7:TestRate=5:TrainingMethod=BFGS" );
// // NN (Multilayer perceptron) - ROOT version
// factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN","!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3" );
// // NN (Multilayer perceptron) - ALEPH version (depreciated)
// factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN","!H:!V:NCycles=2000:HiddenLayers=N+1,N" );
// // Support Vector Machine
// factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001" );
// // Boosted Decision Trees with adaptive boosting
// factory->BookMethod( TMVA::Types::kBDT, "BDT","!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
// // Boosted Decision Trees with gradient boosting
// factory->BookMethod( TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50:NNodesMax=5" );
// // Boosted Decision Trees with bagging
// factory->BookMethod( TMVA::Types::kBDT, "BDTB","!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
// // Predictive learning via rule ensembles (RuleFit)
// factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit","H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
  factory->TrainAllMethods(); 
  factory->TestAllMethods(); 
  factory->EvaluateAllMethods(); 

  tmvaOut.Write(); 

}

// divdeOf20 by default will run all 20 jobs. But you can just run one job by define it.
void doMLP(int data_start = 123, int data_stop=152, int divideOf20 = 20){
  // doMLP_pair(data_start,data_stop,0,34,0);
  // return;
  int index = 0;
  for(int pol = 0; pol < 2; pol++){
    for(int ant1 = 0; ant1 < 47; ant1++){
      for(int ant2 = ant1+1; ant2 < 48; ant2++){
        if((ant2-ant1+16+2)%16 <= 4){
          //there are 770 pairs for h and v. So divide the work into 20 parts, which can be run paralle.
          if(divideOf20!=20 && index%20 != divideOf20){
            index++;
            continue;
          }

          std::cout<<pol << " " <<ant1<< " "<< ant2<< " \n"<< index<<std::endl;
          // doMLP_pair(data_start,data_stop,ant1,ant2,pol);
          index++;
        }
      }
    }
  }
}


