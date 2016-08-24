#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "BasicFilters.h" 
#include "FilterStrategy.h"
#include "Util.h"
#include "TTree.h"
#include "TFile.h"
#include "UCFilters.h"
#include "AnalysisConfig.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"


void doSimulated(int run = 1, int max = 0, bool sine_subtract = false)
{

  FFTtools::loadWisdom("wisdom.dat"); 

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN; 

  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 2; 
  

  UCorrelator::Analyzer analyzer(&cfg); 

  TString outname; 
  if (max) outname.Form("simulation/simulated_%d_max_%d%s.root",run,max, sine_subtract ? "_sinsub" : "" ); 
  else outname.Form("simulation/simulated_%d%s.root",run, sine_subtract ? "_sinsub" : "" ); 

  TFile ofile(outname, "RECREATE"); 
  TTree * tree = new TTree("simulation"," Simulated events"); 
  AnitaEventSummary * sum = new AnitaEventSummary; 

  
  AnitaGeomTool *geomTool = AnitaGeomTool::Instance();
  
  double sourceLat, sourceLon, sourceAlt, sourceMag;
  double thetaWave,phiWave;
  int inu;
  TString icemcfilename;
  icemcfilename.Form("$ANITA_ROOT_DATA/run%d/icefinal%d.root", run, run);
  TFile infile(icemcfilename, "READ"); 
  TTree *icetree = (TTree*)infile.Get("passing_events");
  icetree->SetBranchAddress("inu",          &inu              );
  icetree->SetBranchAddress("sourceLon",    &sourceLon        );
  icetree->SetBranchAddress("sourceLat",    &sourceLat        );
  icetree->SetBranchAddress("sourceMag",    &sourceMag        );
  
 
  FilterStrategy strategy (&ofile); 
  if (sine_subtract) 
  {
    double fmins[1] = {0.2}; 
    double fmaxs[1] = {1.3}; 
    strategy.addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4,1,fmins,fmaxs)); 
    strategy.addOperation(new SimplePassBandFilter(0.2,1.3)); 
    strategy.addOperation(new ALFAFilter); 
  }

  printf("Strategy applied!\n"); 

  RawAnitaHeader *hdr = 0 ; 
  Adu5Pat *patptr = 0; 
  tree->Branch("summary",          &sum      ); 
  tree->Branch("header",           &hdr      ); 
  tree->Branch("pat",              &patptr   );
  tree->Branch("thetaExpectedDeg", &thetaWave);
  tree->Branch("phiExpectedDeg",   &phiWave  );

  int ndone = 0; 


  for (int i =0 ; i < d.N(); i++) {

    d.getEntry(i); 
    printf("----(%d)-----\n",i);
    
    UsefulAdu5Pat pat(d.gps()); 
    
    printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

    analyzer.analyze(&ev, sum); 
    ofile.cd(); 
    hdr = d.header(); 
    patptr = d.gps(); 

    icetree->GetEntry(i);

    if (hdr->eventNumber!=inu){
      std::cout << " We have a problem with eventNumbers : " << hdr->eventNumber << " " << inu << std::endl;
      break;
    }
    sourceAlt=sourceMag-geomTool->getDistanceToCentreOfEarth(sourceLat);

    pat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaWave,phiWave);
    thetaWave*=TMath::RadToDeg();
    phiWave*=TMath::RadToDeg();
    
    std::cout << " Theta wave IceTree: " << thetaWave << std::endl;
    std::cout << " Phi wave IceTree: " << phiWave << std::endl;



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
   
  int run = nargs < 2 ? 1 : atoi(args[1]); 
  int max = nargs < 3 ? 0 : atoi(args[2]); 
  int sinsub = nargs < 4 ? 0 : atoi(args[3]); 

  doSimulated(run,max,sinsub); 


}
