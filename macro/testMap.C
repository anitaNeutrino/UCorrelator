#include "ProbabilityMap.h" 
#include "UCFilters.h" 
#include "PointingResolutionModel.h"
#include "Analyzer.h" 
#include "AnalysisConfig.h" 
#include "AnitaDataset.h" 
#include "SystemResponse.h" 


void testMap(int run =342, int event = 58023120, AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal, int ith = 0) 
{

  FFTtools::loadWisdom("wisdom.dat"); 
  AnitaDataset d(run); 
  d.getEvent(event); 


  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg); 
  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("sinsub_5_3_ad_3"), d.gps(), d.header()); 
  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum,d.truth()); 

  StereographicGrid g(1024,1024); 
  UCorrelator::ConstantPointingResolutionModel m(0.2,0.3);
  UCorrelator::ProbabilityMap map(&g,&m); 

  TFile f("debugfile.root","RECREATE"); 
  map.add(&sum,d.gps(), pol, ith,1,&f); 

  map.segmentationScheme()->Draw("colz", map.getProbabilities()); 
  double igl = 0; 

  for (int i = 0; i < g.NSegments(); i++)
    igl += map.getProbabilities()[i]; 

  printf("integral: %g\n", igl); 

  FFTtools::saveWisdom("wisdom.dat"); 
}
