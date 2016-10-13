#include "Analyzer.h" 
#include "AnitaDataset.h"
#include "BasicFilters.h"
#include "AnalysisConfig.h"
#include "UCFilters.h"
#include "SystemResponse.h"
#include "FFTtools.h"

UCorrelator::Analyzer *doInteractive(int run = 352, int event = 60849734, bool decimated = false )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy strategy; 
 
  TF1 fn("foo"," (x < 0.2) * exp((x-0.2)/0.01)  + (x > 0.2 && x < 1.2) * (1-0.05*x) + (x > 1.2) * exp((1.2-x)/0.02)", 0,2); 


  strategy.addOperation(new UCorrelator::SineSubtractFilter);
//  strategy.addOperation(new SimplePassBandFilter(0.2,1.2)); 
  strategy.addOperation(new ALFAFilter); 

  AnitaDataset d(run,decimated);

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 

  FilteredAnitaEvent ev(d.useful(),&strategy, d.gps(), d.header()); 



  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter; 
  cfg.deconvolution_method = new UCorrelator::WienerDeconvolution(&fn); 
//  cfg.response_option = UCorrelator::AnalysisConfig::ResponseHarmSignalOnly; 
  //cfg.combine_unfiltered = false; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 

  AnitaEventSummary sum; 
  analyzer->analyze(&ev,&sum); 
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


  FFTtools::saveWisdom("wisdom.dat"); 
  return analyzer; 
}
