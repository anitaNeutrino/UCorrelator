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


void doldb(int run = 352, int max = 0, bool sine_subtract = false)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max) outname.Form("ldb/%d_max_%d%s.root",run,max, sine_subtract ? "_sinsub" : "" ); 
  else outname.Form("ldb/%d%s.root",run, sine_subtract ? "_sinsub" : "" ); 

  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("ldb","ldb Hpol"); 
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

    if (UCorrelator::isLDBHPol(&pat, d.header()) || UCorrelator::isLDBVPol(&pat,d.header()))
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

int main (int nargs, char ** args)
{
   
  int run = nargs < 2 ? 352 : atoi(args[1]); 
  int max = nargs < 3 ? 0 : atoi(args[2]); 
  int sinsub = nargs < 4 ? 0 : atoi(args[3]); 

  doldb(run,max,sinsub); 


}
