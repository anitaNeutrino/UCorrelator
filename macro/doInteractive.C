#include "Analyzer.h" 
#include "AnitaDataset.h"
#include "BasicFilters.h"
#include "AnalysisConfig.h"
#include "UCFilters.h"
#include "FFTtools.h"

UCorrelator::Analyzer *doInteractive(int run = 352, int event = 60832108, bool decimated = false )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy strategy; 
 // strategy.addOperation(new SimplePassBandFilter(0.2,1.3)); 
 
  strategy.addOperation(new UCorrelator::SineSubtractFilter);
  strategy.addOperation(new ALFAFilter); 

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
