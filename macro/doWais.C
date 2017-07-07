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

void doWais(int run = 352, int max = 0, const char * filter = "sinsub_10_3_ad_2")
{

  FFTtools::loadWisdom("wisdom.dat"); 


  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  cfg.start_pol = AnitaPol::kHorizontal; 
  cfg.end_pol = AnitaPol::kHorizontal; 
  

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max && start) outname.Form("wais/%d_max_%d_start_%d_%s.root",run,max,start,filter); 
  else if (max) outname.Form("wais/%d_max_%d_%s.root",run,max,filter); 
  else if (start) outname.Form("wais/%d_start_%d_%s.root",run,start,filter); 
  else outname.Form("wais/%d_%s.root",run, filter); 


  TFile ofile(outname, "RECREATE"); 

  TTree * tree = new TTree("wais","WAIS Hpol"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 
  UCorrelator::fillStrategyWithKey(&strategy, filter); 

  RawAnitaHeader *hdr = 0 ; 
  UsefulAdu5Pat *patptr = 0; 
  tree->Branch("summary",&sum); 
  tree->Branch("header",&hdr); 
  tree->Branch("pat",&patptr); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
  {

    d.getEntry(i); 
    printf("----(%d)-----\n",i); 

    UsefulAdu5Pat pat(d.gps()); 

    if (UCorrelator::isWAISHPol(&pat, d.header()))
    {
      printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
      FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

      analyzer.analyze(&ev, sum); 
      ofile.cd(); 
      hdr = d.header();
      patptr = &pat; 
      tree->Fill(); 
      ndone++; 
    }

    if (max && ndone > max) break; 

  }

  ofile.cd(); 
  tree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}
