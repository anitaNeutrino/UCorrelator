<<<<<<< HEAD
#include "AnitaDataset.h"
#include "FFTtools.h" 
#include "AnalysisConfig.h"
=======

#include "FFTtools.h" 
>>>>>>> changespecavg

UCorrelator::Analyzer *doInteractive(int run = 342, int event = 58023120, bool decimated = false, bool simulated = false )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy *strategy= new FilterStrategy; 
 
//  TF1 fn("foo"," (x < 0.2) * exp((x-0.2)/0.01)  + (x > 0.2 && x < 1.2) * (1-0.05*x) + (x > 1.2) * exp((1.2-x)/0.02)", 0,2); 

//  UCorrelator::SpectrumAverage *avg = new UCorrelator::SpectrumAverage(run,60); 
//  avg->computePeakiness(); 

 // UCorrelator::CombinedSineSubtractFilter * ssf = new UCorrelator::CombinedSineSubtractFilter(0.05,10);
 // ssf->setInteractive(true); 
//  printf("UCorrelator::CombinedSineSubtractFilter* ssf = (UCorrelator::CombinedSineSubtractFilter*) %p\n", ssf); 
//  strategy.addOperation(ssf);
//  strategy.addOperation(new SimplePassBandFilter(0.2,1.2)); 

// UCorrelator::SineSubtractFilter * ssf = new UCorrelator::SineSubtractFilter(0.10,3);
//  ssf->makeAdaptive(avg); 
//  UCorrelator::AdaptiveMinimumPhaseFilter * mp = new UCorrelator::AdaptiveMinimumPhaseFilter(avg,-2,5); 
//  printf("UCorrelator::AdaptiveMinimumPhaseFilter * mp = (UCorrelator::AdaptiveMinimumPhaseFilter *) %p\n",mp); 
//  strategy->addOperation(mp); 
//
//UCorrelator::AdaptiveButterworthFilter * butter = new UCorrelator::AdaptiveButterworthFilter(&avg); 
//printf("UCorrelator::AdaptiveButterworthFilter * butter = (UCorrelator::AdaptiveButterworthFilter *) %p\n",butter); 
//strategy.addOperation(butter); 


  AnitaDataset d(run,decimated, WaveCalType::kDefault, simulated ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA, AnitaDataset::kRandomizePolarity );

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 



  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 2; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  //cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  cfg.enable_group_delay= !simulated; 
  cfg.delay_to_center = true; 
//  cfg.r_time_shift_correction = !simulated; 
//  cfg.max_peak_trigger_angle =45; 
//  cfg.correlator_theta_lowest=90; 
//  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
//  cfg.response_option = UCorrelator::AnalysisConfig::ResponseHarmSignalOnly; 
  //cfg.combine_unfiltered = false; 
 

  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 
  cfg.use_coherent_spectra = false; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 

//  strategy->addOperation(new UCorrelator::DeconvolveFilter(analyzer->getResponseManager(), cfg.deconvolution_method)); 
//  strategy->addOperation(new ALFAFilter); 
//  strategy->addOperation(ssf); 

//  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("adsinsub_3_10_3"), d.gps(), d.header()); 
  UCorrelator::setAdaptiveFilterSpectrumAverageNSecs(10); 

  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"), d.gps(), d.header()); 

  printf("auto fae = (FilteredAnitaEvent *) %p;\n",ev); 

  //ev->plotSummary(0,0); 


  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum,d.truth()); 
  TCanvas * chpol = new TCanvas("chpol","hpol",1800,1000); 
  TCanvas * cvpol = new TCanvas("cvpol","vpol",1800,1000); 
  analyzer->drawSummary( chpol,cvpol, true); 
  /*
  TCanvas * c2 = new TCanvas; 
  c2->Divide(2,1); 
  const UCorrelator::AbstractResponse * response = analyzer->getResponseManager()->response(0,0); 
  AnalysisWaveform * imp =response->impulseResponse(0.1, 981); 
  c2->cd(1); 
  imp->drawEven(); 
  AnalysisWaveform * deconv = response->deconvolve(imp,cfg.deconvolution_method); 
  c2->cd(2); 
  deconv->drawEven(); 
  */

  /*
  TCanvas * filter =  new TCanvas("filter","Filter"); 
  filter->Divide(2,1); 
  filter->cd(1); 
  mp->getCurrentFilterPower(AnitaPol::kHorizontal,0)->Draw(); 
  filter->cd(2); 
  mp->getCurrentFilterTimeDomain(AnitaPol::kHorizontal,0)->Draw(); 
  */

   

//  butter->getFilter(AnitaPol::kHorizontal,0)->drawResponse(0,101,10); 

  FFTtools::saveWisdom("wisdom.dat"); 


  AnitaTemplateSummary ats; 
  AnitaTemplateMachine *atm = new AnitaTemplateMachine; 
  atm->loadTemplates(); 

  atm->doTemplateAnalysis( analyzer->getCoherent(AnitaPol::kHorizontal,0),0,0,&ats);
  atm->doTemplateAnalysis( analyzer->getCoherent(AnitaPol::kVertical,0),1,0,&ats);

  printf("cray4[0][0]: %g\n", ats.coherent[0][0].cRay[4]); 
  printf("cray4[1][0]: %g\n", ats.coherent[1][0].cRay[4]); 


  return analyzer; 
}
