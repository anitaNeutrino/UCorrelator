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
  cfg.enable_group_delay = false;
  

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
  double thetaWave2,phiWave2;
  double weight;
  double posnu[3];
  double rfexit[5][3];
  double r_bn[3], r_enterice[3], r_in[3];
  int inu;
  TString icemcfilename;
  icemcfilename.Form("$ANITA_ROOT_DATA/run%d/icefinal%d.root", run, run);
  TFile infile(icemcfilename, "READ"); 
  TTree *icetree = (TTree*)infile.Get("passing_events");
  icetree->SetBranchAddress("inu",          &inu              );
  icetree->SetBranchAddress("sourceLon",    &sourceLon        );
  icetree->SetBranchAddress("sourceLat",    &sourceLat        );
  icetree->SetBranchAddress("sourceAlt",    &sourceAlt        );
  icetree->SetBranchAddress("sourceMag",    &sourceMag        );
  icetree->SetBranchAddress("posnu",        &posnu            );
  icetree->SetBranchAddress("rfexit",       &rfexit           );
  icetree->SetBranchAddress("r_bn",         &r_bn             );
  icetree->SetBranchAddress("r_enterice",   &r_enterice       );
  icetree->SetBranchAddress("r_in",         &r_in             );
  icetree->SetBranchAddress("weight",       &weight           );

  
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
  tree->Branch("summary",           &sum       ); 
  tree->Branch("header",            &hdr       ); 
  tree->Branch("pat",               &patptr    );
  tree->Branch("thetaExpectedDeg",  &thetaWave );
  tree->Branch("phiExpectedDeg",    &phiWave   );
  tree->Branch("thetaExpectedDeg2", &thetaWave2);
  tree->Branch("phiExpectedDeg2",   &phiWave2  );
  tree->Branch("weight",            &weight    );

  int ndone = 0; 
  double tempLon, tempLat, tempAlt;
  
  for (int i =0 ; i < d.N(); i++) {
  // for (int i =0 ; i < 1; i++) {

    d.getEntry(i); 
    printf("----(%d)-----\n",i);
    
    UsefulAdu5Pat pat(d.gps()); 
    
    printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 

    analyzer.analyze(&ev, sum); 
    ofile.cd(); 
    hdr = d.header(); 
    patptr = d.gps(); 

    std::cout << " Measured phi and theta : " << sum->peak[1][0].phi << " " << sum->peak[1][0].theta << std::endl;

    icetree->GetEntry(i);

    if (hdr->eventNumber!=inu){
      std::cout << " We have a problem with eventNumbers : " << hdr->eventNumber << " " << inu << std::endl;
      break;
    }

    pat.getThetaAndPhiWave(sourceLon, sourceLat, sourceAlt, thetaWave,phiWave);
    thetaWave*=TMath::RadToDeg();
    phiWave*=TMath::RadToDeg();
    
    std::cout << " Theta wave IceTree: " << thetaWave << std::endl;
    std::cout << " Phi wave IceTree: " << phiWave << std::endl;
    
    geomTool->getLatLonAltFromCartesian(posnu, tempLat, tempLon, tempAlt);
    pat.getThetaAndPhiWave(tempLon, tempLat, tempAlt, thetaWave2,phiWave2);
    thetaWave2*=TMath::RadToDeg();
    phiWave2*=TMath::RadToDeg();
    std::cout << " Theta wave posnu: " << thetaWave2 << std::endl;
    std::cout << " Phi wave posnu: " << phiWave2 << std::endl;

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
