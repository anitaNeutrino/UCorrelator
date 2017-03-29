
#include "macro/cuts.C"
const int IWAIS=0;
const int IND_WAIS =0;








void doPulserPlots(const char * pulser, TChain *c)
{

  const char * filter = c->GetName(); 

  TCanvas * hpol = new TCanvas(TString::Format("%s_hpol_%s", pulser,filter),TString::Format("%s H-Pol with %s",pulser,filter), 1800,900); 
  hpol->Divide(2,1); 
  hpol->cd(1); 
  c->Draw(TString::Format("FFTtools::wrap(peak[0][0].theta-%s.theta,360,0):FFTtools::wrap(peak[0][0].phi-%s.phi,360,0)",pulser,pulser), brightestPeak && "coherent[0][0].totalPower > 2*coherent[0][0].totalPowerXpol"); 
  hpol->cd(2); 
  TH2 * h_hpol = new TH2D(TString::Format("h_hpol_%s_%s", pulser, filter ),"Zoomed H-Pol; dphi; dtheta;",100,-3,3,100,-3,3); 

  c->Draw(TString::Format("FFTtools::wrap(peak[0][0].theta-%s.theta,360,0):FFTtools::wrap(peak[0][0].phi-%s.phi,360,0) >> h_hpol_%s_%s",pulser,pulser,pulser,filter), brightestPeak && "coherent[0][0].totalPower > 2*coherent[0][0].totalPowerXpol","colz"); 

  hpol->SaveAs(TString::Format("filterPlots/%s_hpol_%s.pdf",pulser,filter)); 

  TCanvas * vpol = new TCanvas(TString::Format("%s_vpol_%s", pulser,filter),TString::Format("%s V-Pol with %s",pulser,filter), 1800,900); 
  vpol->Divide(2,1); 
  vpol->cd(1); 
  c->Draw(TString::Format("FFTtools::wrap(peak[1][0].theta-%s.theta,360,0):FFTtools::wrap(peak[1][0].phi-%s.phi,360,0)",pulser,pulser), brightestPeak && "coherent[1][0].totalPower > 1.25*coherent[1][0].totalPowerXpol"); 
  vpol->cd(2); 
  TH2 * h_vpol = new TH2D(TString::Format("h_vpol_%s_%s", pulser, filter ),"Zoomed V-Pol; dphi; dtheta",100,-3,3,100,-3,3); 
  c->Draw(TString::Format("FFTtools::wrap(peak[1][0].theta-%s.theta,360,0):FFTtools::wrap(peak[1][0].phi-%s.phi,360,0) >> h_vpol_%s_%s",pulser,pulser,pulser,filter), brightestPeak && "coherent[1][0].totalPower > 1.25*coherent[1][0].totalPowerXpol ","colz"); 
  vpol->SaveAs(TString::Format("filterPlots/%s_vpol_%s.pdf",pulser,filter)); 


  TCanvas * dpol = new TCanvas(TString::Format("%s_dpol_%s", pulser,filter),TString::Format("%s D-Pol with %s",pulser,filter), 1800,900); 
  dpol->Divide(2,1); 
  dpol->cd(1); 
  c->Draw(TString::Format("FFTtools::wrap(peak[][0].theta-%s.theta,360,0):FFTtools::wrap(peak[][0].phi-%s.phi,360,0)",pulser,pulser), brightestPeak && "abs(coherent[][0].totalPower-coherent[][0].totalPowerXpol)/coherent[][0].totalPower  < 0.25"); 
  dpol->cd(2); 
  TH2 * h_dpol = new TH2D(TString::Format("h_dpol_%s_%s", pulser, filter ),"Zoomed D-Pol; dphi; dtheta",100,-3,3,100,-3,3); 
  c->Draw(TString::Format("FFTtools::wrap(peak[][0].theta-%s.theta,360,0):FFTtools::wrap(peak[][0].phi-%s.phi,360,0) >> h_dpol_%s_%s",pulser,pulser,pulser,filter), brightestPeak && "abs(coherent[][0].totalPower-coherent[][0].totalPowerXpol)/coherent[][0].totalPower  < 0.25","colz"); 
  dpol->SaveAs(TString::Format("filterPlots/%s_dpol_%s.pdf",pulser,filter)); 

  

  TCanvas * fraction_filtered = new TCanvas(TString::Format("%s_%s_filter_fraction",pulser,filter), TString::Format("%s Filter Fraction with %s", pulser,filter),1800,900); 

  fraction_filtered->Divide(2,2); 
  fraction_filtered->cd(1); 

  TH2 * h_fracangle = new TH2D(TString::Format("h_fracangle_%s_%s", pulser, filter ),"Filtered Fraction by Angle (all); dphi ; dtheta",100,-10,10,100,-10,10); 
  c->Draw(TString::Format("FFTtools::wrap(peak[][].theta - %s.theta,360,0): FFTtools::wrap(peak[][].phi-%s.phi,360,0):coherent_filtered[][].totalPower / coherent[][].totalPower  >> h_fracangle_%s_%s",pulser,pulser,pulser,filter),brightestPeak,"colz"); 
  TH2 * h_fractime = new TH2D(TString::Format("h_fractime_%s_%s", pulser, filter ),"Filtered Fraction by Angle (all); triggerTime ; triggerTimeNs",100,c->GetMinimum("triggerTime"),c->GetMaximum("triggerTime"),1000,0,1e9); 
  h_fractime->GetXaxis()->SetTimeDisplay(1); 
  fraction_filtered->cd(2); 
  c->Draw(TString::Format("triggerTimeNs:triggerTime:coherent_filtered[][].totalPower / coherent[][].totalPower  >> h_fractime_%s_%s",pulser,filter),brightestPeak,"colz"); 

  fraction_filtered->cd(3); 
  c->Draw(TString::Format("coherent_filtered[][].totalPower / coherent[][].totalPower: coherent.totalPower"),brightestPeak,"colz"); 

  fraction_filtered->cd(4); 
  c->Draw(TString::Format("coherent_filtered[][].totalPower / coherent[][].totalPower: coherent.peakHilbert"),brightestPeak,"colz"); 

  fraction_filtered->SaveAs(TString::Format("filterPlots/%s_fraction_%s.pdf",pulser,filter)); 

  TCanvas * distance_plot = new TCanvas(TString::Format("%s_%s_distance",filter,pulser),TString::Format("Distance Plot for %s %s",filter,pulser)); 

  TH2 * h_distance = new TH2D(TString::Format("h_distance_%s_%s", pulser, filter ),"Angular distance (all); triggerTime ; triggerTimeNs",1000,c->GetMinimum("triggerTime"),c->GetMaximum("triggerTime"),1000,0,1e9); 
  c->Draw(TString::Format("triggerTimeNs:triggerTime:sqrt(pow(FFTtools::wrap(peak[][].theta - %s.theta,360,0),2) + pow(FFTtools::wrap(peak[][].phi-%s.phi,360,0),2))  >> h_distance_%s_%s",pulser,pulser,pulser,filter),brightestPeak,"colz"); 
  distance_plot->SaveAs(TString::Format("filterPlots/%s_distance_%s.pdf",pulser,filter)); 
}



void makeCutPlot(const char * filter, TChain * bg, TChain * wais, TChain * ldb) 
{

  gStyle->SetOptStat(0); 
  TCanvas * c = new TCanvas(TString::Format("ccut_%s",filter), TString::Format("Cutplot %s",filter)); 
  TH2I * hbg = new TH2I(TString::Format("hbg_%s", filter),TString::Format("The standard cut plot for %s; Correlation Map Peak; Coherent Peak Hilbert",filter), 100,0,0.5,100,0,200); 
  bg->Draw(TString::Format("coherent[][].peakHilbert:peak[][].value >> hbg_%s",filter),brightestPeak && aboveHorizon && blastCut,"colz"); 

  wais->SetMarkerColor(2); 
  wais->Draw("coherent[][].peakHilbert:peak[][].value", brightestPeak && "abs(FFTtools::wrap(peak[][].phi-wais.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[][].theta - wais.theta,3600)) < 5","psame"); 

  ldb->SetMarkerColor(3); 
  ldb->Draw("coherent[][].peakHilbert:peak[][].value", brightestPeak && "abs(FFTtools::wrap(peak[][].phi-ldb.phi,360,0)) < 5 && abs(FFTtools::wrap(peak[][].theta - ldb.theta,3600)) < 5","psame"); 

  c->SaveAs(TString::Format("filterPlots/standard_cut_plot_%s.pdf",filter)); 
  gStyle->SetOptStat(1); 

}




void doFilterAlgo(const char * filter) 
{
  TChain cwais(filter); 
  cwais.Add("filter/*_wais_*.root"); 

  doPulserPlots("wais",&cwais); 


  TChain cldb(filter); 
  cldb.Add("filter/*_ldb_*.root"); 

  doPulserPlots("ldb",&cldb); 

  TChain cbg(filter); 
  cbg.Add("filter/*_bg*.root"); 

  makeCutPlot(filter, &cbg,&cwais,&cldb); 

}

void filterEvalA3Plots()
{
  system("mkdir -p filterPlots"); 

  doFilterAlgo("sinsub"); 
  doFilterAlgo("adsinsub"); 
  doFilterAlgo("butter"); 
  doFilterAlgo("minphase"); 

}
