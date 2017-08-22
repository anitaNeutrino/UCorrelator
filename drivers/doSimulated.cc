#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "TruthAnitaEvent.h" 
#include "BasicFilters.h" 
#include "SystemResponse.h" 
#include "FilterStrategy.h"
#include "UCUtil.h"
#include "TTree.h"
#include "TFile.h"
#include "UCFilters.h"
#include "AnalysisConfig.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"


void doSimulated(int run = 1, int max = 0, int start = 0, const char * filter = "sinsub_10_3_ad_2")
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,false,WaveCalType::kDefault,AnitaDataset::ANITA_MC_DATA); // Monte Carlo! 
  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 3;
  cfg.enable_group_delay = false; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::ImpulseResponseXCorr; 
  cfg.max_peak_trigger_angle = 90; 

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 

  if (max && start) outname.Form("simulated/%d_max_%d_start_%d_%s.root",run,max,start,filter); 
  else if (max) outname.Form("simulated/%d_max_%d_%s.root",run,max,filter); 
  else if (start) outname.Form("simulated/%d_start_%d_%s.root",run,start,filter); 
  else outname.Form("simulated/%d_%s.root",run, filter); 




  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("simulation"," Simulated events"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 

  
  FilterStrategy strategy (&ofile); 
  UCorrelator::fillStrategyWithKey(&strategy, filter); 
  printf("Strategy applied!\n"); 

  RawAnitaHeader *hdr = 0 ; 
  Adu5Pat *patptr = 0; 
  TruthAnitaEvent * truth = 0; 
  tree->Branch("summary",           &sum       ); 
  tree->Branch("header",            &hdr       ); 
  tree->Branch("pat",               &patptr    );
  tree->Branch("truth",               &truth    );
  int ndone = 0; 
  
  for (int i =start ; i < d.N(); i++) {
  // for (int i =0 ; i < 1; i++) {

    d.getEntry(i); 
    printf("----(%d)-----\n",i);
    
    UsefulAdu5Pat pat(d.gps()); 
    
    printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

    analyzer.analyze(&ev, sum,d.truth()); 

    ofile.cd(); 
    hdr = d.header(); 
    patptr = d.gps(); 
    truth = d.truth(); 

    tree->Fill(); 

    ndone++; 

    if (max && ndone >= max) break; 

  }

  ofile.cd(); 
  tree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}

int main (int nargs, char ** args)
{
   
  int run = nargs < 2 ? 352 : atoi(args[1]); 
  int max = nargs < 3 ? 0 : atoi(args[2]); 
  int start = nargs < 4 ? 0 : atoi(args[3]); 
  const char * filter = nargs < 5 ? 0 :args[4]; 


  if (filter) 
    doSimulated(run,max,start,filter); 
  else
    doSimulated(run,max,start); 


}
