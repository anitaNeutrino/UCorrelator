#include "ProbabilityMap.h" 


void testMap(int run =352, int event = 5802310, AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal, int ith = 0) 
{

  AnitaDataset(d,false); 
  d.getEvent(event); 

  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new UCorrelator::AllPassDeconvolution; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg); 
  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("sinsub_5_3_ad_3"), d.gps(), d.header()); 
  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum,d.truth()); 

  ProbabilityMap map; 
  map.add(&sum,d.gps(), pol, ith); 

  map->segmentationScheme()->Draw("colz", map.getProbabilities()); 

}
