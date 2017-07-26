#include "FFTtools.h"
#include "Analyzer.h"
#include "TF1.h" 
#include "FilteredAnitaEvent.h"
#include "SystemResponse.h" 
#include "FilterStrategy.h"
#include "UCUtil.h"
#include "TTree.h"
#include "TFile.h"
#include "UCFilters.h"
#include "BasicFilters.h" 
#include "AnalysisConfig.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"

#ifdef UCORRELATOR_OPENMP
#include <omp.h>
#endif 

void doAll(int run = 352, int max = 0, int start = 0, const char * filter = "sinsub_10_3_ad_2")
{

  FFTtools::loadWisdom("wisdom.dat"); 

  AnitaVersion::set(3); 

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  
#ifdef UCORRELATOR_OPENMP
  printf("Max threads: %d\n", omp_get_max_threads()); 
#endif


  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution(); 

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max && start) outname.Form("a%dall/%d_max_%d_start_%d_%s.root",AnitaVersion::get(),run,max,start,filter); 
  else if (max) outname.Form("a%dall/%d_max_%d_%s.root",AnitaVersion::get(),run,max,filter); 
  else if (start) outname.Form("a%dall/%d_start_%d_%s.root",AnitaVersion::get(),run,start,filter); 
  else outname.Form("a%dall/%d_%s.root",AnitaVersion::get(),run, filter); 


  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree(TString::Format("anita%d",AnitaVersion::get()),TString::Format("anita%d", AnitaVersion::get())); 
  tree->SetAutoFlush(1000); 
  AnitaEventSummary * sum = new AnitaEventSummary; 

  FilterStrategy strategy (&ofile); 

  UCorrelator::fillStrategyWithKey(&strategy, filter); 

//  printf("Strategy applied!\n"); 

  RawAnitaHeader *hdr = 0 ; 
  Adu5Pat *patptr = 0; 
  tree->Branch("summary",&sum); 
  tree->Branch("header",&hdr); 
  tree->Branch("pat",&patptr); 

  int ndone = 0; 
  for (int i =start ; i < d.N(); i++)
  {

    d.getEntry(i); 
    printf("----(%d)-----\n",i); 


    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

    analyzer.analyze(&ev, sum); 
    ofile.cd(); 
    hdr = d.header(); 
    patptr = d.gps(); 
    tree->Fill(); 
    ndone++; 

    if (max && ndone > max) break; 

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
    doAll(run,max,start,filter); 
  else
    doAll(run,max,start); 

  

}
