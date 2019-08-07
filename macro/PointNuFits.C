
#ifndef CUTS_LOADED
#include "macro/cuts.C"
#endif

#include "AnitaTMVA.h" 



void makePointingTrees()
{
  TChain c("simulation"); 
  c.Add("simulated_hooper/*.root"); 

  TFile of("pointing_training.root","RECREATE"); 

  AnitaTMVA::MVAVarSet phi_vars("phi_point.tmva"); 
  AnitaTMVA::MVAVarSet theta_vars("theta_point.tmva"); 


  TTree * dphi = makeTMVATree(&c, &of, "dphi",phi_vars,""); 
  TTree * dtheta = makeTMVATree(&c, &of, "dtheta",theta_vars,""); 

  dphi->Write();
  dtheta->Write(); 
}




void doPointing() 
{
  TChain dphic("dphi"); 
  TChain dthetac("dtheta"); 
  dphic.Add("pointing_training.root");
  dthetac.Add("pointing_training.root");

  TFile fout("pointingTMVA.root","RECREATE"); 

  TMVA::Factory * f= new TMVA::Factory("pointing",&fout,"V:Color:DrawProgressBar:AnalysisType=Regression"); 

  AnitaTMVA::MVAVarSet theta_vars("theta_point.tmva"); 
  AnitaTMVA::MVAVarSet phi_vars("phi_point.tmva"); 

  TMVA::DataLoader * dl_phi = new TMVA::DataLoader("dphi_training"); 
  TMVA::DataLoader * dl_theta = new TMVA::DataLoader("dtheta_training"); 

  dl_phi->AddRegressionTree(&dphic);
  dl_theta->AddRegressionTree(&dthetac);

  dl_phi->SetWeightExpression("weight","Regression");
  dl_theta->SetWeightExpression("weight","Regression");

  theta_vars.setUpData(dl_theta);
  phi_vars.setUpData(dl_phi);

  TCut cut("(isMostImpulsive && pointsToMC && weight > 1e-6)"); 
  dl_phi->AddCut(cut); 
  dl_theta->AddCut(cut); 

  f->BookMethod(dl_phi, TMVA::Types::kLD, "LD",""); 
  f->BookMethod(dl_theta, TMVA::Types::kLD, "LD",""); 
  f->BookMethod(dl_phi, TMVA::Types::kKNN, "KNN","nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim"); 
  f->BookMethod(dl_theta, TMVA::Types::kKNN, "kNN","nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim"); 
  f->BookMethod(dl_phi, TMVA::Types::kMLP, "MLP",""); 
  f->BookMethod(dl_theta, TMVA::Types::kMLP, "MLP",""); 


  f->TrainAllMethods(); 
  f->TestAllMethods(); 
  f->EvaluateAllMethods(); 
  fout.Write(); 

}
