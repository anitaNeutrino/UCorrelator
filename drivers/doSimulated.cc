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
#include "BH13Filter.h"
#include "Hical2.h"



void doSimulated(int run = 1, int max = 0, int start = 0, const char * out_dir = "simulated", const char * filter = "sinsub_10_3_ad_2")
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,false,WaveCalType::kDefault,AnitaDataset::ANITA_MC_DATA); // Monte Carlo! 
  UCorrelator::AnalysisConfig cfg; 
    cfg.nmaxima = 3;
    cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;
    cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 

  if (max && start) outname.Form("%s/%d_max_%d_start_%d_%s.root",out_dir,run,max,start,filter); 
  else if (max) outname.Form("%s/%d_max_%d_%s.root",out_dir,run,max,filter); 
  else if (start) outname.Form("%s/%d_start_%d_%s.root",out_dir,run,start,filter); 
  else outname.Form("%s/%d_%s.root",out_dir,run, filter); 


  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("simulation"," Simulated events"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 

  double dtheta = 1.; double dphi = 2.; bool blockout = true;
  analyzer.setTrackSun(dtheta, dphi, blockout);

  FilterStrategy* forDeco = new FilterStrategy;
  forDeco->addOperation(new UCorrelator::AntiBH13Filter());
  analyzer.setExtraFiltersDeconvolved(forDeco);
  analyzer.setDisallowedAntennas(0, (1ul<<45));  // Vpol ant45 is bad! So disable it.

  FilterStrategy strategy (&ofile);
  UCorrelator::fillStrategyWithKey(&strategy, filter);
  strategy.addOperation(new UCorrelator::BH13Filter());

  RawAnitaHeader *hdr = 0 ; 
  Adu5Pat *patptr = 0; 
  double isHC = 0;
  // TruthAnitaEvent * truth = 0; 
  tree->Branch("summary",           &sum       ); 
  tree->Branch("header",            &hdr       ); 
  tree->Branch("pat",               &patptr    );
  tree->Branch("isHC",              &isHC    );
  // tree->Branch("truth",               &truth    );
  int ndone = 0; 
  
  for (int i =start ; i < d.N(); i++) {
  // for (int i =0 ; i < 1; i++) {

    d.getEntry(i); 
    printf("----(%d)-----\n",i);
    //fixed run46 hkfile problem.
    // if(d.header()->realTime >= 1480725529 and d.header()->realTime<=1480730678){
    //   printf("Skip this event which from run 46.\n");
    //   continue;
    // }
    
    UsefulAdu5Pat pat(d.gps()); 
    const time_t ctt = time(0);
    printf("Processing event %d (%d) \t|%s", d.header()->eventNumber,ndone,asctime(localtime(&ctt)));
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

    analyzer.analyze(&ev, sum,d.truth()); 

    ofile.cd(); 
    hdr = d.header(); 
    patptr = d.gps();
    isHC = Hical2::isHical(sum); 
    // truth = d.truth(); 

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
  const char * outdir = nargs < 5 ? "simulated" : args[4]; 
  const char * filter = nargs < 6 ? 0 :args[5]; 


  if (filter) 
    doSimulated(run,max,start,outdir, filter); 
  else
    doSimulated(run,max,start,outdir); 


}
