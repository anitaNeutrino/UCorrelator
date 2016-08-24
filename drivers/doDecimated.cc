#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "FilterStrategy.h"
#include "Util.h"
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

void doDecimated(int run = 352, int max = 0, bool sine_subtract = false)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,true); 
  UCorrelator::AnalysisConfig cfg; 
  
#ifdef UCORRELATOR_OPENMP
  printf("Max threads: %d\n", omp_get_max_threads()); 
#endif

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max) outname.Form("decimated/%d_max_%d%s.root",run,max, sine_subtract ? "_sinsub" : "" ); 
  else outname.Form("decimated/%d%s.root",run, sine_subtract ? "_sinsub" : "" ); 

  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("decimated","Decimated"); 
  tree->SetAutoFlush(1000); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 
  if (sine_subtract) 
  {
    double fmins[1] = {0.2}; 
    double fmaxs[1] = {1.3}; 
    strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4,1,fmins,fmaxs)); 
    strategy.addOperation(new SimplePassBandFilter(0.2,1.3)); 
    strategy.addOperation(new ALFAFilter); 
  }
  else
  {
    UCorrelator::applyAbbysFilterStrategy(&strategy); 
  }

//  printf("Strategy applied!\n"); 

  RawAnitaHeader *hdr = 0 ; 
  Adu5Pat *patptr = 0; 
  tree->Branch("summary",&sum); 
  tree->Branch("header",&hdr); 
  tree->Branch("pat",&patptr); 

  int ndone = 0; 
  for (int i =0 ; i < d.N(); i++)
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
  int sinsub = nargs < 4 ? 1 : atoi(args[3]); 

  doDecimated(run,max,sinsub); 


}
