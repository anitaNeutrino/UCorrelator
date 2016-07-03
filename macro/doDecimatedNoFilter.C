#include "FFTtools.h" 
#include "AnitaDataset.h" 
#include "Analyzer.h" 

void doDecimatedNoFilter(int run = 352, int max = 0)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,true); 

  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 5; 
  
  UCorrelator::Analyzer* analyzer = new UCorrelator::Analyzer(&cfg); 

  TFile ofile(TString::Format("decimated_no_filter/%d.root", run), "RECREATE"); 
  TTree * tree = new TTree("decimated_no_filter","decimated_no_filter"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 
  FilterStrategy strategy; 
  strategy.addOperation(new UCorrelator::SimplePassBandFilter(0.2,1.3)); 


  tree->Branch("summary",&sum); 
  RawAnitaHeader *hdr = 0;
  Adu5Pat *pat =  0; 
  tree->Branch("header",&hdr); 
  tree->Branch("pat",&pat); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
  {

  printf("----(%d)-----\n",i); 
    d.getEntry(i); 
    printf("%d\n",i); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 
    analyzer->analyze(&ev, sum); 
    ofile.cd(); 
    header = d.header(); 
    pat = d.gps(); 
    tree->Fill(); 
    ndone++; 

    if (max && ndone > max) break; 

  }

  ofile.cd(); 
  tree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}
