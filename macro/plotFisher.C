void plotFisher(const char * var = "F")
{
  TChain csim("simulation"); 
  csim.Add("thermalTrees/simulated*.root"); 

  TChain cbg("anita3"); 
  cbg.Add("thermalTrees/a3all*.root"); 
  cbg.SetProof(1); 

  csim.SetLineColor(4); 
  cbg.SetLineColor(2); 

  TFile fisher_plots("fisher_plots.root","RECREATE");
  TH1D * hbg = new TH1D("hbg","Upward-Pointing Thermals That Pass Quality Cuts",200,-10,10); 
  TH1D * hsim = new TH1D("hsim","Weighted MC Pass Quality Cuts",200,-10,10); 

  TCanvas * dists = new TCanvas(); 

  cbg.Draw("F >> hbg",
              "isMostImpulsive * (chisq < 0.01) * (MaxPeak < 1000) * (theta < 0)  * (!payloadBlast) *  ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))", 
             ""); 

  csim.Draw("F >> hsim",
             "weight * isMostImpulsive * (chisq < 0.01) * (MaxPeak < 1000) * pointsToMC * (!payloadBlast) *  ((iteration < 5 && HPolTrigger) || (iteration > 4 && VPolTrigger))", 
             "same"); 


  fisher_plots.cd(); 
  hbg->Write(); 
  hsim->Write(); 


}
