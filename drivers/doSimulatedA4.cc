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

// Set run to run number used in icemc
// Requires new A4 spec averages, see the script in ../data
// Run without _ad_2 in filter to avoid using default terminated thermal spectrum for A3 for peakiness
void doSimulatedA4(int run = 4, int max = 0, int start = 0, const char * out_dir = "simulated", const char * filter = "sinsub_10_3")
{
  cout << "Processing simulated data for ANITA-4" << endl;

  FFTtools::loadWisdom("wisdom.dat");
  AnitaVersion::set(4);

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run,false,WaveCalType::kDefault,AnitaDataset::ANITA_MC_DATA); // Monte Carlo! 
  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 3;
  cfg.enable_group_delay = false;
  // Added option
  cfg.use_forced_trigger_rms = false; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
//  cfg.max_peak_trigger_angle = 90; 
  cfg.fill_blast_fraction = true; 
  cfg.combine_nantennas = 15; 
  cfg.zoomed_nant = 15; 

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname;

  if (max && start) outname.Form("%s/%d_max_%d_start_%d_%s.root",out_dir,run,max,start,filter); 
  else if (max) outname.Form("%s/%d_max_%d_%s.root",out_dir,run,max,filter); 
  else if (start) outname.Form("%s/%d_start_%d_%s.root",out_dir,run,start,filter); 
  else outname.Form("%s/%d_%s.root",out_dir,run, filter); 

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
  double percDone = 0;
  
  for (int i =start ; i < d.N(); i++) {
  // for (int i =0 ; i < 1; i++) {

    percDone = double(i)/double(d.N()) * 100;
    
    d.getEntry(i); 
    
    UsefulAdu5Pat pat(d.gps()); 
    
    printf("Processing event %d",d.header()->eventNumber);
    printf("-----(%d)-----(%f %%)\n",i,percDone);
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
  const char * outdir = nargs < 5 ? "simulated" : args[4]; 
  const char * filter = nargs < 6 ? 0 :args[5]; 


  if (filter) 
    doSimulatedA4(run,max,start,outdir, filter); 
  else
    doSimulatedA4(run,max,start,outdir); 


}
