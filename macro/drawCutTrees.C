void drawCutTrees(const char * filter="sinsub_10_3_ad_2")
{

  TFile *f = new TFile (TString::Format("thermalCutTrees_%s.root",filter)); 

  TTree * signal = (TTree*) f->Get("signal_in"); 
  TTree * bg = (TTree*) f->Get("bg_in"); 
  signal->SetMarkerColorAlpha(3,0.05); 
  signal->SetMarkerStyle(1); 

  TCanvas * c = new TCanvas(filter,filter,1000,1000); 
 
  c->Divide(2,2); 

  c->cd(1); 
  gPad->SetLogz(); 
  bg->Draw("coherentHilbertPeak:mapPeak >> h1(100,0,0.5,100,0,300)","","colz"); 
  signal->Draw("coherentHilbertPeak:mapPeak","","psame"); 

  c->cd(2); 
  bg->Draw("coherentHilbertPeak:mapPeak >> h2(100,0,0.5,100,0,300)","dPhiSun > 10","colz"); 
  signal->Draw("coherentHilbertPeak:mapPeak","","psame"); 
  gPad->SetLogz(); 

  c->cd(3); 
  bg->Draw("coherentHilbertPeak:mapPeak >> h3(100,0,0.5,100,0,300)","dPhiNorth > 90","colz"); 
  signal->Draw("coherentHilbertPeak:mapPeak","","psame"); 
  gPad->SetLogz(); 

  c->cd(4); 
  bg->Draw("coherentHilbertPeak:mapPeak >> h4(100,0,0.5,100,0,300)","dPhiNorth > 90 && dPhiSun > 10","colz"); 
  signal->Draw("coherentHilbertPeak:mapPeak","","psame"); 
  gPad->SetLogz(); 

  c = new TCanvas("deconv","deconv",1000,1000); 
 
  c->Divide(2,2); 

  c->cd(1); 
  gPad->SetLogz(); 
  bg->Draw("deconvHilbertPeak:mapPeak >> h1(100,0,0.5,100,0,300)","","colz"); 
  signal->Draw("deconvHilbertPeak:mapPeak","","psame"); 

  c->cd(2); 
  bg->Draw("deconvHilbertPeak:mapPeak >> h2(100,0,0.5,100,0,300)","dPhiSun > 10","colz"); 
  signal->Draw("deconvHilbertPeak:mapPeak","","psame"); 
  gPad->SetLogz(); 

  c->cd(3); 
  bg->Draw("deconvHilbertPeak:mapPeak >> h3(100,0,0.5,100,0,300)","dPhiNorth > 90","colz"); 
  signal->Draw("deconvHilbertPeak:mapPeak","","psame"); 
  gPad->SetLogz(); 

  c->cd(4); 
  bg->Draw("deconvHilbertPeak:mapPeak >> h4(100,0,0.5,100,0,300)","dPhiNorth > 90 && dPhiSun > 10","colz"); 
  signal->Draw("deconvHilbertPeak:mapPeak","","psame"); 
  gPad->SetLogz(); 


}
