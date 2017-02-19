#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "AnitaDataset.h" 
#include "RawAnitaHeader.h"
#include "AnalysisConfig.h"
#include "UCFilters.h"
#include "BasicFilters.h"
#include "Util.h"
#include "FilterStrategy.h"

void doDecimated(int run = 352, int max = 0, bool sine_subtract = false)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,true); 
  UCorrelator::AnalysisConfig cfg; 
  

  UCorrelator::Analyzer analyzer(&cfg); 

  TFile ofile(TString::Format("decimated/%d.root", run), "RECREATE"); 
  TTree * tree = new TTree("decimated","decimated"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 

  if (sine_subtract) 
  {
    double fmins[1] = {0.2}; 
    double fmaxs[1] = {1.3}; 
    strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4,1,fmins,fmaxs)); 
    strategy.addOperation(new SimplePassBandFilter(0.2, 1.3)); 
//    strategy.addOperation(new UCorrelator::SimplePassBandFilter(0.2.1.3)); 
  }
  else
  {
    UCorrelator::applyAbbysFilterStrategy(&strategy); 
  }


//  printf("Strategy applied!\n"); 

  tree->Branch("summary",&sum); 
  RawAnitaHeader *hdr = 0; 
  Adu5Pat *pat = 0; 
//  RawAnitaHeader *hdr; 
//  Adu5Pat *pat; 
  tree->Branch("header",&hdr); 
  tree->Branch("pat",&pat); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
  {

    printf("----(%d)-----\n",i); 
    d.getEntry(i); 
//    printf("%d\n",i); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 
    analyzer.analyze(&ev, sum); 
    ofile.cd(); 
    hdr = d.header();
//    header = d.header(); 
    pat = d.gps(); 
    tree->Fill(); 
    ndone++; 

    if (max && ndone > max) break; 

  }

  ofile.cd(); 
  tree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}
