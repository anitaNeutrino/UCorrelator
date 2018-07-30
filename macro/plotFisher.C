// root plotFisher.C
// The output will be in fisher_plot.root

void fillEfficeincy(TH1D * h, TH1D *eff,double scale, const char * cut = "", const char * ordering = ">>") 
{
  for (int i = 1; i <=h->GetNbinsX(); i++){
    eff->SetBinContent(i, h->Integral(ordering[0] == '>' ? i : 0, ordering[0] == '>' ? h->GetNbinsX()+1 : i) / scale ); 
  }

  eff->GetXaxis()->SetTitle("impulsivity");
  eff->GetYaxis()->SetTitle("efficiency"); 
  eff->Write();
  return;
}
void fillLogEfficeincy(TH1D * h, TH1D *eff,double scale, const char * cut = "", const char * ordering = ">>") 
{
  for (int i = 1; i <=h->GetNbinsX(); i++){
    double efficiency_val = h->Integral(ordering[0] == '>' ? i : 0, ordering[0] == '>' ? h->GetNbinsX()+1 : i) / scale;
    if (efficiency_val == 0){
      continue;
    }
    std::cout << i << " "<<efficiency_val << std::endl;
    // eff->SetBinContent(i, log10(efficiency_val) - expect_eff );
    // double expect_eff = log10(3.555e15) - log10(2.731)*71.86 * eff->GetXaxis()->GetBinCenter(i);

    eff->SetBinContent(i, log10(efficiency_val)); 
  }
  // double factor = [0]*exp(-[1]*x
  eff->GetXaxis()->SetTitle("impulsivity");
  eff->GetYaxis()->SetTitle("log10(efficiency)"); 
  eff->Write();
  return;
}

void fillEfficeincyDiff(TH1D * h, TH1D *eff,double scale, const char * cut = "", const char * ordering = ">>") 
{
  for (int i = 1; i <=h->GetNbinsX(); i++){
    double efficiency_val = h->Integral(ordering[0] == '>' ? i : 0, ordering[0] == '>' ? h->GetNbinsX()+1 : i) / scale;
    double imp = h->GetXaxis()->GetBinCenter(i);
    if (efficiency_val == 0 or  imp < 0.65 or imp > 0.74){
      continue;
    }
    // double expect_eff = -26.74*imp + 12.56;
    double expect_eff = log10(3.555e15) - log10(2.731)*71.86 *imp;
    std::cout << log10(efficiency_val) - expect_eff << std::endl;
    eff->Fill(log10(efficiency_val) - expect_eff ); 
  }
  // double factor = [0]*exp(-[1]*x
  eff->GetXaxis()->SetTitle("log(efficiency) - log(fitted_efficiency)");
  eff->GetYaxis()->SetTitle("count"); 
  eff->Write();
  return;
}
void plotFisher(const char * var = "deconvImpulsivity")
{
  TChain chain1("anita4"); 
  chain1.Add("thermalTrees/a4all*10000003*.root");
  TChain chain2("anita4"); 
  chain2.Add("thermalTrees/a4all*10000003*.root");
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

  // 2d histogram plot
  TH2D * hist2Dh1 = new TH2D("hist2Dh1","Above Horizontal Events ",200,0,1.1,200,0,1.1); 
  TH2D * hist2Dh2 = new TH2D("hist2Dh2","Below Horizon Events  ",200,0,1.1,200,0,1.1); 
  TH2D * hist2Dh3 = new TH2D("hist2Dh3","MinBias Energy222 MC",200,0,1.1,200,0,1.1); 
  TH2D * hist2Dh4 = new TH2D("hist2Dh4","Wais data",200,0,1.1,200,0,1.1); 
  chain1.Draw("impulsivityV:impulsivityH >> hist2Dh1","theta>0", "colz"); 
  chain2.Draw("impulsivityV:impulsivityH >> hist2Dh2","theta<=-5.8", "colz"); 
  chain3.Draw("impulsivityV:impulsivityH >> hist2Dh3", "weight", "colz"); 
  chain4.Draw("impulsivityV:impulsivityH >> hist2Dh4","","colz");
  hist2Dh1->GetXaxis()->SetTitle("impulsivityH");
  hist2Dh2->GetXaxis()->SetTitle("impulsivityH");
  hist2Dh3->GetXaxis()->SetTitle("impulsivityH");
  hist2Dh4->GetXaxis()->SetTitle("impulsivityH");
  hist2Dh1->GetYaxis()->SetTitle("impulsivityV");
  hist2Dh2->GetYaxis()->SetTitle("impulsivityV");
  hist2Dh3->GetYaxis()->SetTitle("impulsivityV");
  hist2Dh4->GetYaxis()->SetTitle("impulsivityV");

// efficiency plot
  TH1D * eff1 = new TH1D("eff1","Above Horizontal Events ",200,0,1.1); 
  TH1D * eff2 = new TH1D("eff2","Below Horizon Events  ",200,0,1.1); 
  TH1D * eff3 = new TH1D("eff3","MinBias Energy222 MC",200,0,1.1); 
  TH1D * eff4 = new TH1D("eff4","Wais data",200,0,1.1); 
  TH1D * logEff1 = new TH1D("logEff1","Above Horizontal Events ",200,0,1.1);
  TH1D * logRelativeEff1 = new TH1D("logRelativeEff1","Above Horizontal Events ",50,-1,1);
  double scale1 = h1->Integral();
  double scale2 = h2->Integral();
  double scale3 = h3->Integral();
  double scale4 = h4->Integral();
  fillEfficeincy(h1,eff1,scale1);
  fillEfficeincy(h2,eff2,scale2);
  fillEfficeincy(h3,eff3,scale3);
  fillEfficeincy(h4,eff4,scale4);
  fillLogEfficeincy(h1,logEff1,scale1);
  fillEfficeincyDiff(h1,logRelativeEff1 , scale1);

  fisher_plots.cd(); 
  h1->Write(); 
  h2->Write(); 
  h3->Write(); 
  h4->Write();

  hist2Dh1->Write(); 
  hist2Dh2->Write(); 
  hist2Dh3->Write(); 
  hist2Dh4->Write();

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
  auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
  // legend2->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend2->AddEntry("eff1","Above Horizontal Events","l");
  legend2->AddEntry("eff2","Below horizon Events","l");
  legend2->AddEntry("eff3","MinBias Energy222 MC","l");
  legend2->AddEntry("eff4","Wais data","l");
  legend2->Draw();


  // fit exponetial line

// // exponential fit to impulsivity
gStyle->SetOptFit(0101);
TF1 *expo1 = new TF1("expo1","[0]*exp(-[1]*x)",0.5,0.68);
expo1->SetParNames("A","B");
expo1->SetParameters(3.0e15,70);
eff1->Fit("expo1","","",0.65, 0.8)

// // gaussian fit to impulsivity histogram
gStyle->SetOptFit(0101);
TF1 *gaussianFit = new TF1("gaussianFit","[0]*exp(-[1]*(x-[2])^2 )",0.4,0.8);
gaussianFit->SetParNames("A","B","C");
gaussianFit->SetParameters(2e5,200, 0.5);
h1->Fit("gaussianFit","","",0.45, 0.75)

// // linear fit
// gStyle->SetOptFit(1111);
TF1 *linearFit3 = new TF1("linearFit3","[0]*x + [1]",0.5,0.68);
linearFit3->SetParNames("A","B");
linearFit3->SetParameters(-10,20);
logEff1->Fit("linearFit3","","",0.66, 0.73)
}



