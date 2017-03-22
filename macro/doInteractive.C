#include "Analyzer.h" 
#include "AnitaDataset.h"
#include "BasicFilters.h"
#include "AnalysisConfig.h"
#include "UCFilters.h"
#include "SystemResponse.h"
#include "FFTtools.h"
#include "DigitalFilter.h" 

UCorrelator::Analyzer *doInteractive(int run = 342, int event = 58023120, bool decimated = false )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy strategy; 
 
  TF1 fn("foo"," (x < 0.2) * exp((x-0.2)/0.01)  + (x > 0.2 && x < 1.2) * (1-0.05*x) + (x > 1.2) * exp((1.2-x)/0.02)", 0,2); 

  UCorrelator::SpectrumAverage avg(run,60); 
  avg.computePeakiness(); 

 // UCorrelator::CombinedSineSubtractFilter * ssf = new UCorrelator::CombinedSineSubtractFilter(0.05,10);
 // ssf->setInteractive(true); 
//  printf("UCorrelator::CombinedSineSubtractFilter* ssf = (UCorrelator::CombinedSineSubtractFilter*) %p\n", ssf); 
//  strategy.addOperation(ssf);
//  strategy.addOperation(new SimplePassBandFilter(0.2,1.2)); 

//  UCorrelator::SineSubtractFilter * ssf = new UCorrelator::SineSubtractFilter(0.05,3);
//  ssf->makeAdaptive(&avg); 
//  strategy.addOperation(ssf); 
  UCorrelator::AdaptiveMinimumPhaseFilter * mp = new UCorrelator::AdaptiveMinimumPhaseFilter(&avg,-2,5); 
  printf("UCorrelator::AdaptiveMinimumPhaseFilter * mp = (UCorrelator::AdaptiveMinimumPhaseFilter *) %p\n",mp); 
  strategy.addOperation(mp); 
//
//UCorrelator::AdaptiveButterworthFilter * butter = new UCorrelator::AdaptiveButterworthFilter(&avg); 
//printf("UCorrelator::AdaptiveButterworthFilter * butter = (UCorrelator::AdaptiveButterworthFilter *) %p\n",butter); 
//strategy.addOperation(butter); 

  strategy.addOperation(new ALFAFilter); 

  AnitaDataset d(run,decimated);

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 

  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),&strategy, d.gps(), d.header()); 

  ev->plotSummary(0,0); 



  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter; 
  cfg.deconvolution_method = new UCorrelator::WienerDeconvolution(&fn); 
//  cfg.response_option = UCorrelator::AnalysisConfig::ResponseHarmSignalOnly; 
  //cfg.combine_unfiltered = false; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 

  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum); 
  analyzer->drawSummary(); 
  TCanvas * c2 = new TCanvas; 
  c2->Divide(3,1); 
  const UCorrelator::AbstractResponse * response = analyzer->getResponseManager()->response(0,0); 
  AnalysisWaveform * imp =response->impulseResponse(0.1, 981); 
  c2->cd(1); 
  imp->drawEven(); 
  AnalysisWaveform * deconv = response->deconvolve(imp,cfg.deconvolution_method); 
  c2->cd(2); 
  deconv->drawEven(); 
  c2->cd(3); 
  fn.DrawCopy(); 

  TCanvas * filter =  new TCanvas("filter","Filter"); 
  filter->Divide(2,1); 
  filter->cd(1); 
  mp->getCurrentFilterPower(AnitaPol::kHorizontal,0)->Draw(); 
  filter->cd(2); 
  mp->getCurrentFilterTimeDomain(AnitaPol::kHorizontal,0)->Draw(); 

   

//  butter->getFilter(AnitaPol::kHorizontal,0)->drawResponse(0,101,10); 

  FFTtools::saveWisdom("wisdom.dat"); 
  return analyzer; 
}
