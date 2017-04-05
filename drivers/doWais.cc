#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "BasicFilters.h" 
#include "TF1.h" 
#include "FilterStrategy.h"
#include "SystemResponse.h" 
#include "UCUtil.h"
#include "TTree.h"
#include "TFile.h"
#include "UCFilters.h"
#include "AnalysisConfig.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"


void doWais(int run = 352, int max = 0, bool deconvolve = true)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  cfg.start_pol = AnitaPol::kHorizontal; 
  cfg.end_pol = AnitaPol::kHorizontal; 
  
  if (deconvolve)
  {
    TF1 *fn = new TF1("foo"," (x < 0.2) * exp((x-0.2)/0.01)  + (x > 0.2 && x < 1.2) * (1-0.05*x) + (x > 1.2) * exp((1.2-x)/0.02)", 0,2); 
    cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter; 
    cfg.deconvolution_method = new UCorrelator::WienerDeconvolution(fn); 
  }


  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max) outname.Form("wais/wais_hpol_%d_max_%d%s.root",run,max, deconvolve ? "_deconv" : "" ); 
  else outname.Form("wais/wais_hpol_%d%s.root",run, deconvolve ? "_deconv" : "" ); 

  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("wais","WAIS Hpol"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 


  FilterStrategy strategy (&ofile); 
  double fmins[1] = {0.2}; 
  double fmaxs[1] = {1.3}; 
  strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 1,fmins,fmaxs)); 
  strategy.addOperation(new ALFALanczosFilter); 

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

int main (int nargs, char ** args)
{
   
  int run = nargs < 2 ? 352 : atoi(args[1]); 
  int max = nargs < 3 ? 0 : atoi(args[2]); 
  int sinsub = nargs < 4 ? 1 : atoi(args[3]); 

  doWais(run,max,sinsub); 


}
