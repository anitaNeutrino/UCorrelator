
#include "FFTtools.h" 
#include "AnitaConventions.h" 



int prettyPlot(int event = 15717147, int anitaVer = 3,
               int pol = 0, int pk = 0, 
               bool simulated = false, bool filtered = true, bool true_deconvolve = true, bool freq_instead =true,
               bool blind_polarity = false, bool plot_zoomed = false) 
{

  AnitaVersion::set(anitaVer);  

  FFTtools::loadWisdom("wisdom.dat"); 
  int run = AnitaDataset::getRunContainingEventNumber(event); 
 
  AnitaDataset d(run,false, WaveCalType::kDefault,
                 simulated ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA, 
                 blind_polarity ? AnitaDataset::kRandomizePolarity : AnitaDataset::kDefault);


  event > 0 ? d.getEvent(event) : d.getEntry(-event); 

  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1+pk; 

  TF1 snr("snr", " 100 * ( (x < 0.18) * TMath::Exp((x-0.18)/0.002) + (x >= 0.18 && x <= 0.5) + (x > 0.5) * TMath::Exp((0.5-x)/0.01) )", 0,1.3); 

  if (anitaVer == 3) cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter;
  else if (anitaVer == 4) cfg.response_option == UCorrelator::AnalysisConfig::ResponseTUFF;
  cfg.deconvolution_method =  true_deconvolve ? (AnitaResponse::DeconvolutionMethod*) new AnitaResponse::WienerDeconvolution(&snr): (AnitaResponse::DeconvolutionMethod*) new AnitaResponse::AllPassDeconvolution; 

  cfg.enable_group_delay= !simulated; 
  cfg.delay_to_center = true; 
  cfg.r_time_shift_correction = false;//!simulated; 

  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 
  cfg.use_coherent_spectra = false; 

  cfg.correlation_gain_correction = 40; 
  cfg.correlator_theta_lowest = 50; 
  cfg.correlator_theta_highest = 30; 
  cfg.correlator_nphi = 360; 
  cfg.correlator_ntheta = 160; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 


  UCorrelator::setAdaptiveFilterSpectrumAverageNSecs(10); 

  FilteredAnitaEvent ev(d.useful(),UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"), d.gps(), d.header()); 

  AnitaEventSummary sum; 
  analyzer->analyze(&ev,&sum,d.truth()); 
 

  TCanvas * c = new TCanvas("waveforms","Waveforms", 1800,900); 

  AnalysisWaveform * coherent =  new AnalysisWaveform(*analyzer->getCoherent(AnitaPol::AnitaPol_t(pol),pk, filtered)); 
  AnalysisWaveform * deconv = new AnalysisWaveform(*analyzer->getDeconvolved(AnitaPol::AnitaPol_t(pol),pk, filtered)); 
  AnalysisWaveform * coherent_xpol = new AnalysisWaveform(*analyzer->getCoherent(AnitaPol::AnitaPol_t(1-pol),pk, filtered)); 
  AnalysisWaveform * deconv_xpol = new AnalysisWaveform(*analyzer->getDeconvolved(AnitaPol::AnitaPol_t(1-pol),pk, filtered)); 

  /** Taper the edges */ 


//  FFTtools::TukeyWindow tukey(0.98); 
//
//  tukey.apply(deconv->Neven(), deconv->updateEven()->GetY());
//  tukey.apply(deconv_xpol->Neven(), deconv_xpol->updateEven()->GetY());
//
//  deconv->updateEven()->zeroMean(); 
//  deconv_xpol->updateEven()->zeroMean(); 

  c->SetBorderSize(0); 
  c->Divide(2,1) ; 

  c->cd(1)->SetGridy(); 

  TLegend * leg = new TLegend(0.7,0.7,0.9,0.9); 

  if (freq_instead) 
  {

    TGraph * g_freq = new TGraph(*deconv->powerdB()); 
    TGraph * g_freq_x = new TGraph(*deconv_xpol->powerdB()); 


    g_freq->GetYaxis()->SetTitle("Power (arb)"); 
    g_freq->GetXaxis()->SetTitle("GHz"); 
    g_freq->GetXaxis()->SetRangeUser(0,1.2); 
    g_freq->GetYaxis()->SetTitleOffset(1.2); 
    g_freq->SetTitle(true_deconvolve ? "Deconvolved Spectrum" : "Average Spectrum"); 
    g_freq->SetLineColor(pol == 0 ? 38 : 44); 

    g_freq_x->SetLineColor(pol == 0 ? 44 : 38); 
    g_freq->Draw("al"); 
    g_freq_x->Draw("lsame"); 
 
    leg->AddEntry(g_freq,pol == 0 ? "HPol" : "VPol","l"); 
    leg->AddEntry(g_freq_x,pol == 1 ? "HPol" : "VPol","l"); 
  }

  else
  {
    coherent->padFreq(4); 
    coherent_xpol->padFreq(4); 

    TGraph* g_coh = new TGraph(*coherent->even()); 
    g_coh->SetLineWidth(2); 
    g_coh->GetYaxis()->SetTitle("mV"); 
    g_coh->GetXaxis()->SetTitle("ns"); 
    g_coh->GetXaxis()->SetRangeUser(0,100); 
    g_coh->GetYaxis()->SetTitleOffset(1.2); 
    g_coh->SetTitle("Coherently-Averaged"); 
    g_coh->SetLineColor(pol == 0 ? 38 : 44); 

    TGraph* g_coh_x = new TGraph(*coherent_xpol->even()); 
    g_coh_x->SetLineColor(pol == 0 ? 44 : 38); 
    leg->AddEntry(g_coh,pol == 0 ? "HPol" : "VPol","l"); 
    leg->AddEntry(g_coh_x,pol == 1 ? "HPol" : "VPol","l"); 
    g_coh->Draw("al"); 
    g_coh_x->Draw("lsame"); 
  }

  leg->Draw(); 

  c->cd(2)->SetGridy(); 

  deconv->padFreq(4); 
  deconv_xpol->padFreq(4); 

  TGraph* g_dec = new TGraph(*deconv->even()); 
  g_dec->SetLineWidth(2); 
  g_dec->GetYaxis()->SetTitle(true_deconvolve ? "E-Field (arb)" : "mV"); 
  g_dec->GetYaxis()->SetTitleOffset(1.2); 
  g_dec->GetXaxis()->SetTitle("ns"); 
  g_dec->GetXaxis()->SetRangeUser(0,100); 
  g_dec->SetLineColor(pol == 0 ? 38 : 44); 
  g_dec->SetTitle( true_deconvolve ? "Coherently-Deconvolved": "Coherently-Dedispersed"); 
  g_dec->Draw("al"); 

  TGraph* g_dec_x = new TGraph(*deconv_xpol->even()); 
  g_dec_x->SetLineColor(pol == 0 ? 44 : 38); 
  g_dec_x->Draw("lsame"); 

  leg->Draw(); 

  TCanvas * c2 = new TCanvas("maps","Maps",1800,900); 

  if (plot_zoomed) 
  {
    c2->Divide(2,1); 
    c2->cd(1); 
  }
  ((UCorrelator::gui::Map*) analyzer->getCorrelationMap(AnitaPol::AnitaPol_t(pol)))->SetTitle(TString::Format("Ev. %d, %s Polarization", event, pol == 0 ? "Horizontal" : "Vertical")); 
  ((UCorrelator::gui::Map*) analyzer->getCorrelationMap(AnitaPol::AnitaPol_t(pol)))->GetXaxis()->SetTitle("Payload #phi (degrees)"); 
  ((UCorrelator::gui::Map*) analyzer->getCorrelationMap(AnitaPol::AnitaPol_t(pol)))->GetYaxis()->SetTitle("Elevation (degrees)"); 
  ((UCorrelator::gui::Map*) analyzer->getCorrelationMap(AnitaPol::AnitaPol_t(pol)))->GetYaxis()->SetTitleOffset(1.2); 
  ((UCorrelator::gui::Map*) analyzer->getCorrelationMap(AnitaPol::AnitaPol_t(pol)))->Draw("nh ns colz"); 
  if (plot_zoomed)
  {
    c2->cd(2); 
    ((UCorrelator::gui::Map*) analyzer->getZoomedCorrelationMap(AnitaPol::AnitaPol_t(pol),pk))->SetTitle("Zoomed Map"); 
    ((UCorrelator::gui::Map*) analyzer->getZoomedCorrelationMap(AnitaPol::AnitaPol_t(pol),pk))->GetXaxis()->SetTitle("Payload #phi (degrees)"); 
    ((UCorrelator::gui::Map*) analyzer->getZoomedCorrelationMap(AnitaPol::AnitaPol_t(pol),pk))->GetYaxis()->SetTitle("Elevation (degrees)"); 
    ((UCorrelator::gui::Map*) analyzer->getZoomedCorrelationMap(AnitaPol::AnitaPol_t(pol),pk))->GetYaxis()->SetTitleOffset(1.2); 
    ((UCorrelator::gui::Map*) analyzer->getZoomedCorrelationMap(AnitaPol::AnitaPol_t(pol),pk))->Draw("nh ns colz"); 
  }

  return 0; 
}
