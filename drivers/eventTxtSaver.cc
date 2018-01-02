#include "AnitaDataset.h" 
#include "Analyzer.h" 
#include "AnitaEventSummary.h" 
#include "AnalysisConfig.h" 
#include <fstream> 
#include "SystemResponse.h" 
#include <ostream> 
#include <istream> 
#include "UCFilters.h" 
#include "TGraph.h" 
#include "UsefulAdu5Pat.h" 
#include "FilteredAnitaEvent.h" 

using namespace std; 

void eventTxtSaver(int run, int eventNumber) {

  //AnitaVersion::set(3);

  //const int run=342; //a good wais run
  AnitaDataset data(run);
  
  UCorrelator::AnalysisConfig *config = new UCorrelator::AnalysisConfig(); 
   //config->response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter
  config->response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  config->deconvolution_method = new AnitaResponse::AllPassDeconvolution;
  config->combine_nantennas = 15;
  config->zoomed_nant=15;

  UCorrelator::Analyzer *analyzer = new UCorrelator::Analyzer(config,true); //interactive needs to be true
  FilterStrategy *fStrat = UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2");
  AnitaEventSummary *eventSummary = new AnitaEventSummary();
  
  
  
  
  TGraph HwaveformC; // coherently summed
  TGraph VwaveformC;
  TGraph HwaveformD; // coherent summed deconvolved
  TGraph VwaveformD; 

  int numSavedWaveforms = 0;
  int entry=0;
  
   cerr << "number of entries= "<< data.N() << endl;
 
    data.getEvent(eventNumber);
    
    UsefulAdu5Pat *usefulGPS = new UsefulAdu5Pat(data.gps());
   // cerr << "entry= "<< entry <<" WAIS bool= "<<UCorrelator::isWAISHPol(usefulGPS,data.header(),config)<< endl;
      FilteredAnitaEvent *filteredEvent = new FilteredAnitaEvent(data.useful(), fStrat, data.gps(), data.header());
      analyzer->analyze(filteredEvent, eventSummary); 
      HwaveformC  =  TGraph(*analyzer->getCoherent(AnitaPol::kHorizontal,0,true)->even());
      VwaveformC  =  TGraph(*analyzer->getCoherentXpol(AnitaPol::kHorizontal,0,true)->even());
      HwaveformD =  TGraph(*analyzer->getDeconvolved(AnitaPol::kHorizontal,0,true)->even());
      VwaveformD =  TGraph(*analyzer->getDeconvolvedXpol(AnitaPol::kHorizontal,0,true)->even());     
      
      analyzer->clearInteractiveMemory();
      delete filteredEvent;

    delete usefulGPS;
   
    char outNameHC[40],outNameVC[40],outNameHD[40],outNameVD[40];
    
    sprintf(outNameHC,"HpolC%d.txt",eventNumber);
    sprintf(outNameVC,"VpolC%d.txt",eventNumber);
    sprintf(outNameHD,"HpolD%d.txt",eventNumber);
    sprintf(outNameVD,"VpolD%d.txt",eventNumber);    

  ofstream outFileHC(outNameHC);
  for (int pt=0; pt<HwaveformC.GetN(); pt++) {
      outFileHC << HwaveformC.GetX()[pt] << "  " << HwaveformC.GetY()[pt] << endl ; 
  }
  outFileHC.close();
  
  ofstream outFileVC(outNameVC);
  for (int pt=0; pt<VwaveformC.GetN(); pt++) {
      outFileVC << VwaveformC.GetX()[pt] << "  " << VwaveformC.GetY()[pt] << endl ; 
  }
  outFileVC.close();
  
  ofstream outFileHD(outNameHD);
  for (int pt=0; pt<HwaveformD.GetN(); pt++) {
      outFileHD << HwaveformD.GetX()[pt] << "  " << HwaveformD.GetY()[pt] << endl ; 
  }
  outFileHD.close();
  
  ofstream outFileVD(outNameVD);
  for (int pt=0; pt<VwaveformD.GetN(); pt++) {
      outFileVD << VwaveformD.GetX()[pt] << "  " << VwaveformD.GetY()[pt] << endl ; 
  }
  outFileVD.close();  

}
/*////--------------------------------*/

int main(int argc, char **argv)
{
   if(argc<2) {
      cout << "Usage: " << endl;
      cout << "\t" << argv[0] << "<run> <event>" << endl << endl;
      return 0;
   }    
   
  
   Int_t run = atoi(argv[1]);
   Int_t event = atoi(argv[2]);
   
   eventTxtSaver(run,event);
   
}
