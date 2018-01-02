void plotFisher(const char * var = "F")
{
  TChain caugsim("simulation"); 
  caugsim.Add("thermalTrees/simulated*.root"); 

  TChain cdecsim("simulation"); 
  cdecsim.Add("thermalTrees/fakenoise*root"); 

  TChain cbg("anita3"); 
  cbg.Add("thermalTrees/a3all*.root"); 
//  cbg.SetProof(1); 

  caugsim.SetLineColor(4); 
  cbg.SetLineColor(2); 
  cdecsim.SetLineColor(3); 

  TFile fisher_plots("fisher_plots.root","RECREATE");
  TH1D * hbg = new TH1D("hbg","Upward-Pointing Thermals That Pass Quality Cuts",200,-10,10); 
  TH1D * haugsim = new TH1D("haugsim","Weighted Aug17 MC Pass Quality Cuts",200,-10,10); 
  TH1D * hdecsim = new TH1D("hdecsim","Weighted Dec17 MC Pass Quality Cuts",200,-10,10); 

  TCanvas * dists = new TCanvas(); 

  cbg.Draw("F >> hbg",
              "isMostImpulsive *  (MaxPeak < 1000) * (theta < 0)  * (!payloadBlast) *  ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))", 
             ""); 

  caugsim.Draw("F >> haugsim",
             "weight * isMostImpulsive * (MaxPeak < 1000) * pointsToMC * (!payloadBlast) *  ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))", 
             "same"); 

  cdecsim.Draw("F >> hdecsim",
             "weight * isMostImpulsive * (MaxPeak < 1000) * pointsToMC * (!payloadBlast) *  ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))", 
             "same"); 



  fisher_plots.cd(); 
  hbg->Write(); 
  hdecsim->Write(); 
  haugsim->Write(); 


}
