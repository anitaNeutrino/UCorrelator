
#include "FFTtools.h" 
#include "AnitaConventions.h" 

int prettyPlot(int event = 58023120,
               int pol = 0, int pk = 0, 
               bool simulated = false, bool filtered = true, bool true_deconvolve = false, 
               bool blind_polarity = false) 
{

  FFTtools::loadWisdom("wisdom.dat"); 
  int run = AnitaDataset::getRunContainingEventNumber(event); 
 
  AnitaDataset d(run,false, WaveCalType::kDefault,
                 simulated ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA, 
                 blind_polarity ? AnitaDataset::kRandomizePolarity : AnitaDataset::kDefault);


  event > 0 ? d.getEvent(event) : d.getEntry(-event); 

  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1+pk; 

  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method =  true_deconvolve ? (AnitaResponse::DeconvolutionMethod*) new AnitaResponse::WienerDeconvolution: (AnitaResponse::DeconvolutionMethod*) new AnitaResponse::AllPassDeconvolution; 

  cfg.enable_group_delay= !simulated; 
  cfg.delay_to_center = true; 
  cfg.r_time_shift_correction = false;//!simulated; 

  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 
  cfg.use_coherent_spectra = false; 

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

  c->Divide(2,1); 

  c->cd(1); 

  coherent->padFreq(4); 
  coherent_xpol->padFreq(4); 

  TGraph* g_coh = new TGraph(*coherent->even()); 
  g_coh->GetYaxis()->SetTitle("mV"); 
  g_coh->GetXaxis()->SetTitle("ns"); 
  g_coh->SetTitle("Coherently-Summed Waveform ( + cross-pol)"); 
  g_coh->Draw("al"); 

  TGraph* g_coh_x = new TGraph(*coherent_xpol->even()); 
  g_coh_x->SetLineColor(11); 
  g_coh_x->Draw("lsame"); 


  c->cd(2); 


  deconv->padFreq(4); 
  deconv_xpol->padFreq(4); 

  TGraph* g_dec = new TGraph(*deconv->even()); 
  g_dec->SetLineWidth(2); 
  g_dec->GetYaxis()->SetTitle(true_deconvolve ? "V/m" : "mV"); 
  g_dec->GetXaxis()->SetTitle("ns"); 
  g_dec->SetTitle( true_deconvolve ? "Coherently-Summed Deconvolved Waveform ( + cross-pol)": "Coherently-Summed Dedispersed Waveform ( + cross-pol)"); 
  g_dec->Draw("al"); 

  TGraph* g_dec_x = new TGraph(*deconv_xpol->even()); 
  g_dec_x->SetLineColor(11); 
  g_dec_x->Draw("lsame"); 


  TCanvas * c2 = new TCanvas("maps","Maps",1800,900); 

  c2->Divide(2,1); 
  c2->cd(1); 
  ((UCorrelator::gui::Map*) analyzer->getCorrelationMap(AnitaPol::AnitaPol_t(pol)))->Draw("nh colz"); 
  c2->cd(2); 
  ((UCorrelator::gui::Map*) analyzer->getZoomedCorrelationMap(AnitaPol::AnitaPol_t(pol),pk))->Draw("nh colz"); 

  return 0; 
}
