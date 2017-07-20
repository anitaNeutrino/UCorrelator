
#include "AnitaEventSummary.h"
#include "FFTtools.h"
#include "AnitaDataset.h" 
#include "AnalysisConfig.h" 
#include "TEllipse.h" 
#include "TMarker.h" 
#include "UCFilters.h" 
#include "TStyle.h" 
#include "PeakFinder.h" 
#include "FilterStrategy.h" 
#include "Correlator.h" 
#include "TFile.h" 
#include "TH2.h" 
#include "Analyzer.h" 
#include "TCanvas.h" 

AnitaEventSummary * testMasking(int run = 125, int entry = 0, bool nofilter = true) 
{

//  AnalysisWaveform::enableDebug(true); 

  FFTtools::loadWisdom("wisdom.dat"); 

  gStyle->SetOptStat(0); 

  AnitaDataset d(run); 
  if (entry > -1)
  {
    d.getEntry(entry); 
  }
  else
  {
    d.getEntry(-entry); 
  }


  TFile out("test.root","RECREATE"); 
  FilterStrategy strategy(&out); 
  double fmins[]={0.23, 0.43}; 
  double fmaxs[]={0.29, 0.49}; 
  
  if (!nofilter) 
  {
    strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4)); 
  }
  else
  {
    strategy.addOperation(new SimplePassBandFilter(0.2,1.3)); 
//    strategy.addOperation(new UCorrelator::SimplePassBandFilter(0.2,1.3)); 
//    UCorrelator::applyAbbysFilterStrategy(&strategy); 
  }

  printf("strategy applied!\n"); 

//  AnalysisWaveform::enableDebug(true); 
  printf("creating event\n"); 
  FilteredAnitaEvent * fae = new FilteredAnitaEvent(d.useful(), &strategy, d.gps(), d.header(),true); 



  printf("processed strategy!\n"); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  UCorrelator::AnalysisConfig cfg; 
  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg, true);
	//analyzer->setExcludePhiThetaRange(100,200,-60,40);
	//analyzer->setTrackSun(10,25);
	analyzer->setTrackWAIS(10,25);
  analyzer->analyze(fae, sum); 

	analyzer->drawSummary();

  FFTtools::saveWisdom("wisdom.dat"); 
  return sum; 

}
