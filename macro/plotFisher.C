// root plotFisher.C
// The output will be in fisher_plot.root
void plotFisher(const char * var = "deconvImpulsivity")
{
  TChain chain1("anita4"); 
  chain1.Add("thermalTrees/a4all*30002*.root");
  TChain chain2("anita4"); 
  chain2.Add("thermalTrees/a4all*30002*.root");
  TChain chain3("simulation"); 
  chain3.Add("thermalTrees/simulated*1001*.root");
  TChain chain4("wais"); 
  chain4.Add("thermalTrees/wais*root"); 

  chain1.SetLineColor(1); 
  chain2.SetLineColor(2); 
  chain3.SetLineColor(3); 
  chain4.SetLineColor(4); 

  TFile fisher_plots("fisher_plots.root","RECREATE");


  // TH1D * h1 = new TH1D("h1","Above Horizontal Thermals ",200,-1,2); 
  // TH1D * h2 = new TH1D("h2","Below Horizontal Thermals",200,-1,2); 
  // TH1D * h3 = new TH1D("h3","Energy222 MC",200,-1,2); 
  // TH1D * h4 = new TH1D("h4","Wais data",200,-1,2);
  // TCanvas * dists = new TCanvas("c1","c1"); 
  // chain1.Draw("F >> h1","theta>0", ""); 
  // chain2.Draw("F >> h2","theta<0", "same"); 
  // chain3.Draw("F >> h3", "weight", "same"); 
  // chain4.Draw("F >> h4","","same"); 

  TH1D * h1 = new TH1D("h1","Above Horizontal Events ",200,0,1.1); 
  TH1D * h2 = new TH1D("h2","Below Horizontal Events  ",200,0,1.1); 
  TH1D * h3 = new TH1D("h3","MinBias Energy222 MC",200,0,1.1); 
  TH1D * h4 = new TH1D("h4","Wais data",200,0,1.1); 
  TCanvas * dists = new TCanvas("c1","c1"); 
  chain1.Draw("deconvImpulsivity >> h1","theta>-3.5", ""); 
  chain2.Draw("deconvImpulsivity >> h2","theta<=-3.5", "same"); 
  chain3.Draw("deconvImpulsivity >> h3", "weight", "same"); 
  chain4.Draw("deconvImpulsivity >> h4","","same");

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
   // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry("h1","theta>-3.5 Events","l");
   legend->AddEntry("h2","theta<-3.5 Events","l");
   legend->AddEntry("h3","MinBias Energy222 MC","l");
   legend->AddEntry("h4","Wais data","l");
   legend->Draw();
   // h1->GetXaxis()->SetTitle("Fisher");
   h1->GetXaxis()->SetTitle("impulsivity");

  fisher_plots.cd(); 
  h1->Write(); 
  h2->Write(); 
  h3->Write(); 
  h4->Write();
  legend->Write();

  std::cout << "finished"<<std::endl;

}
