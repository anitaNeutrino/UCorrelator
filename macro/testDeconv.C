#include "FFTtools.h" 

UCorrelator::Analyzer *testDeconv(int event = 72164985) 
{

  FFTtools::loadWisdom("wisdom.dat"); 
 
  //load a random run.. the dataset will then load the correct run :) 
  AnitaDataset d(200,false); 

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 



  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseCustomString; 
  cfg.response_string = "newA4responses";
  cfg.deconvolution_method = new AnitaResponse::CLEAN; 
  //cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  cfg.enable_group_delay= true; 
  cfg.delay_to_center = true; 
  cfg.r_time_shift_correction = true; 
  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 
  cfg.use_coherent_spectra = true; 
  cfg.cross_correlate_hv = 1; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 
  UCorrelator::setAdaptiveFilterSpectrumAverageNSecs(10); 
  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"), d.gps(), d.header()); 

  // if you want the pointer to the FilteredAnitaEvent, it's printed here 
  printf("auto fae = (FilteredAnitaEvent *) %p;\n",ev); 


  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum,d.truth()); 
  TCanvas * chpol = new TCanvas("chpol","hpol",1800,1000); 
  TCanvas * cvpol = new TCanvas("cvpol","vpol",1800,1000); 
  analyzer->drawSummary( chpol,cvpol, true); 
  FFTtools::saveWisdom("wisdom.dat"); 
  return analyzer; 
}
