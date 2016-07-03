#include "Analyzer.h" 
#include "FFTtools.h"

UCorrelator::Analyzer *doInteractive(int run = 342, int event = 0, bool decimated = true )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy strategy; 
  strategy.addOperation(new UCorrelator::SimplePassBandFilter(0.2,1.3)); 
  AnitaDataset d(run,decimated);

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 

  FilteredAnitaEvent ev(d.useful(),&strategy, d.gps(), d.header()); 

  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 3; 
  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 

  AnitaEventSummary sum; 
  analyzer->analyze(&ev,&sum); 
  analyzer->drawSummary(); 
  FFTtools::saveWisdom("wisdom.dat"); 
  return analyzer; 
}
