#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "AnitaEventCalibrator.h"
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


void doResolutions( int run = 352, int max = 0, int start = 0, const char * filter = "sinsub_10_3_ad_2",  const char * outputDir = "photogrammetry" )
{

  FFTtools::loadWisdom("wisdom.dat"); 

  AnitaVersion::set(3);
  
  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->usePhotogrammetryNumbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(Int_t surf=0; surf<NUM_SURF; surf++){
    for(Int_t chan=0; chan<NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0; ///< From phase center to AMPAs (hopefully)
    }
  }
  
  AnitaPol::AnitaPol_t pol;

  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  Int_t cutTimeNs = 1200;
  char cpol[100];
  bool isLDB = false;

  Int_t delayGenerator = 199996850; // delay generator was not exactly 200 ms
  Int_t delay=0;
  Int_t constdelay = 500;
  Double_t triggerTimeNsExpected;
  Double_t triggerTimeNs;
  Int_t deltaTriggerTimeNs;

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  
  std::string pulser="";
  
  if (run>300){ // WAIS
    pol = AnitaPol::kHorizontal;
    pulser+="WAIS";
    cfg.start_pol = AnitaPol::kHorizontal; 
    cfg.end_pol = AnitaPol::kHorizontal; 
  } else if (run<150){ // LDB VPOL 
    pol = AnitaPol::kVertical;
    isLDB = true;
    delay =  25000000; // V-POL pulse at 25 ms runs 145-149
    pulser+="LDB";
    cfg.start_pol = AnitaPol::kVertical; 
    cfg.end_pol = AnitaPol::kVertical; 
  } else if (run<154){ // LDB HPOL
    pol = AnitaPol::kHorizontal;
    isLDB = true;
    delay =  50000000; // H-POL pulse at 50 ms
    pulser+="LDB";    
    cfg.start_pol = AnitaPol::kHorizontal; 
    cfg.end_pol = AnitaPol::kHorizontal; 
  } else if (run<172){ // LDB VPOL
    pol = AnitaPol::kVertical;
    isLDB = true;
    delay =  50000000; // V-POL pulse at 50 ms runs 154-171
    pulser+="LDB";
    cfg.start_pol = AnitaPol::kVertical; 
    cfg.end_pol = AnitaPol::kVertical; 
  } else {
    std::cout << "Unknown run" << std::endl;
    return;
  }

  if (pol == AnitaPol::kVertical) sprintf(cpol, "VPOL");
  else if (pol == AnitaPol::kHorizontal) sprintf(cpol, "HPOL");
  
  
//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 
  
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution(); 

  
  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  outname.Form("%s/Resolution_%s_%s_%d_%s.root", outputDir, pulser.c_str(), cpol, run, filter); 
  
  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("resolution","resolution"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 

  FilterStrategy strategy (&ofile); 
  UCorrelator::fillStrategyWithKey(&strategy, filter); 

  RawAnitaHeader *hdr = 0 ; 
  UsefulAdu5Pat *patptr = 0; 
  tree->Branch("summary",&sum); 
  tree->Branch("header",&hdr); 
  tree->Branch("pat",&patptr); 

  int ndone = 0;
  
  bool check = false;
  
  for (int i =start ; i < d.N(); i++)
  {

    d.getEntry(i); 
    //    printf("----(%d)-----\n",i); 

    UsefulAdu5Pat pat(d.gps()); 
    
    if (isLDB){
      
      triggerTimeNs         = d.header()->triggerTimeNs; 
      triggerTimeNsExpected = pat.getLDBTriggerTimeNs();
      deltaTriggerTimeNs    = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      deltaTriggerTimeNs    = deltaTriggerTimeNs%(delayGenerator) - delay - constdelay;
      
      check = (TMath::Abs(deltaTriggerTimeNs) < cutTimeNs);
      
    }else{
      check=UCorrelator::isWAISHPol(&pat, d.header());
    }
    
    if (check)
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
  const char * outDir = nargs < 6 ? 0 :args[5];
  
  if (filter) {
    if (outDir)
      doResolutions(run,max,start,filter, outDir);
    else
      doResolutions(run,max,start,filter); 
  }
  else
    doResolutions(run,max,start); 


}
