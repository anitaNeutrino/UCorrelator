// root plotFisher.C
// The output will be in fisher_plot.root

void fillEfficeincy(TH1D * h, TH1D *eff,double scale, const char * cut = "", const char * ordering = ">>") 
{
  for (int i = 1; i <=h->GetNbinsX(); i++){
    eff->SetBinContent(i, h->Integral(ordering[0] == '>' ? i : 0, ordering[0] == '>' ? h->GetNbinsX()+1 : i) / scale ); 
  }

  eff->GetXaxis()->SetTitle("impulsivity");  d
  eff->GetYaxis()->SetTitle("efficiency"); 
  eff->Write();
  return;
}
void plotFisher(const char * var = "deconvImpulsivity")
{
  TChain chain1("anita4"); 
  chain1.Add("thermalTrees/a4all*30002*.root");
  TChain chain2("anita4"); 
  chain2.Add("thermalTrees/a4all*30002*.root");
  TChain chain3("simulation"); 
  chain3.Add("thermalTrees/simulated*1000*.root");
  TChain chain4("wais"); 
  chain4.Add("thermalTrees/wais*root"); 

  chain1.SetLineColor(1); 
  chain2.SetLineColor(2); 
  chain3.SetLineColor(3); 
  chain4.SetLineColor(4); 

  TFile fisher_plots("fisher_plots.root","RECREATE");
// histogram plot
  TH1D * h1 = new TH1D("h1","Above Horizontal Events ",200,0,1.1); 
  TH1D * h2 = new TH1D("h2","Below Horizon Events  ",200,0,1.1); 
  TH1D * h3 = new TH1D("h3","MinBias Energy222 MC",200,0,1.1); 
  TH1D * h4 = new TH1D("h4","Wais data",200,0,1.1); 
  TCanvas * dists = new TCanvas("c1","c1"); 
  chain1.Draw("impulsivity >> h1","theta>0", ""); 
  chain2.Draw("impulsivity >> h2","theta<=-5.8", "same"); 
  chain3.Draw("impulsivity >> h3", "weight", "same"); 
  chain4.Draw("impulsivity >> h4","","same");
  h1->GetXaxis()->SetTitle("impulsivity");
  h2->GetXaxis()->SetTitle("impulsivity");
  h3->GetXaxis()->SetTitle("impulsivity");
  h4->GetXaxis()->SetTitle("impulsivity");
  h1->GetYaxis()->SetTitle("count");
  h2->GetYaxis()->SetTitle("count");
  h3->GetYaxis()->SetTitle("count");
  h4->GetYaxis()->SetTitle("count");

// efficiency plot
  TH1D * eff1 = new TH1D("eff1","Above Horizontal Events ",200,0,1.1); 
  TH1D * eff2 = new TH1D("eff2","Below Horizon Events  ",200,0,1.1); 
  TH1D * eff3 = new TH1D("eff3","MinBias Energy222 MC",200,0,1.1); 
  TH1D * eff4 = new TH1D("eff4","Wais data",200,0,1.1); 
  double scale1 = h1->Integral();
  double scale2 = h2->Integral();
  double scale3 = h3->Integral();
  double scale4 = h4->Integral();
  fillEfficeincy(h1,eff1,scale1);
  fillEfficeincy(h2,eff2,scale2);
  fillEfficeincy(h3,eff3,scale3);
  fillEfficeincy(h4,eff4,scale4);


  fisher_plots.cd(); 
  h1->Write(); 
  h2->Write(); 
  h3->Write(); 
  h4->Write();
  legend->Write();

  std::cout << "finished"<<std::endl;

}

void backup(){
  // legend
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry("h1","Above Horizontal Events","l");
  legend->AddEntry("h2","Below horizon Events","l");
  legend->AddEntry("h3","MinBias Energy222 MC","l");
  legend->AddEntry("h4","Wais data","l");
  legend->Draw();


  // legend
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry("eff1","Above Horizontal Events","l");
  legend->AddEntry("eff2","Below horizon Events","l");
  legend->AddEntry("eff3","MinBias Energy222 MC","l");
  legend->AddEntry("eff4","Wais data","l");
  legend->Draw();

}



