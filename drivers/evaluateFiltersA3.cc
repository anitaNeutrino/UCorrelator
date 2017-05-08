/** 
 *
 * ANITA-3 evaluation framework executable. 
 *
 * Takes a run, classifies it as one of LDB, WAIS or Background (decimated) and runs pointing 
 * with all filters. 
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu>
 *
 */

#include "TruthAnitaEvent.h" 
#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "FilterStrategy.h"
#include "UCUtil.h"
#include "TTree.h"
#include "TFile.h"
#include "UCFilters.h"
#include <vector>
#include <string> 
#include "SystemResponse.h" 
#include "BasicFilters.h" 
#include "GeomFilter.h" 
#include "AnalysisConfig.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"
#include "WaveformCombiner.h" 
#include "SpectrumAverage.h" 



int min_ldb = 130; 
int max_ldb = 164; 

int min_wais = 332; 
int max_wais = 354; 


FilterStrategy justAlfa; 
std::vector<FilterStrategy*> strategies; 
std::vector<std::string> names; 

void addStrategy(const char * key)
{

  strategies.push_back(UCorrelator::getStrategyWithKey(key)); 
  names.push_back(key); 

} 
UCorrelator::Analyzer * analyzer; 
UCorrelator::AnalysisConfig cfg; 


/** add filter strategies here ! */ 

void setupFilters(TFile* out) 
{

  (void) out; 
  justAlfa.addOperation(new ALFAFilter); 


  addStrategy("sinsub_10_0"); 
  addStrategy("adsinsub_1_10_0"); 
  addStrategy("adsinsub_2_10_0"); 
  addStrategy("adsinsub_1_10_3"); 
  addStrategy("adsinsub_2_10_3"); 
  addStrategy("adsinsub_3_10_3"); 
  addStrategy("adsinsub_2_20_0"); 
  addStrategy("brickwall_2_0"); 
  addStrategy("brickwall_2_1"); 
  addStrategy("geom"); 

}




int main(int nargs, char ** args) 
{

  AnitaVersion::set(3); 
  FFTtools::loadWisdom("wisdom.dat"); 
  int run = atoi(args[1]); 

  bool isWAIS = run >=min_wais && run <= max_wais; 
  bool isLDB = run >=min_ldb && run <= max_ldb; 
  bool isMC = run <0; 
  bool isBG  = !isWAIS && !isLDB && !isMC; 

  AnitaDataset d(abs(run), isBG,WaveCalType::kDefault,isMC ? 0 : -1); //only use decimated if background 

  int max = nargs > 2 ? atoi(args[2]) : 0; 

  TString outname; 
  const char * label = isWAIS ? "wais" : isLDB ? "ldb" : isMC ? "mc" :  "bg"; 

  if (max) outname.Form("filter/%d_%s_max_%d.root", run, label, max); 
  else outname.Form("filter/%d_%s_max_%d.root", run, label, max); 

  TFile ofile(outname, "RECREATE"); 

  if (!ofile.IsOpen())
  {
    return 1; 
  }

  cfg.response_option=UCorrelator::AnalysisConfig::ResponseIndividualBRotter;
  cfg.deconvolution_method = new UCorrelator::AllPassDeconvolution; 
  analyzer = new UCorrelator::Analyzer(&cfg); 

  setupFilters(&ofile); 
 
  UCorrelator::WaveformCombiner wc(12,3,true,true,analyzer->getResponseManager()); //this code is so shitty 
  AnitaEventSummary * sum = new AnitaEventSummary; 
  RawAnitaHeader *hdr = 0 ; 
  UsefulAdu5Pat *patptr = 0; 

  std::vector<TTree*> trees(strategies.size()); 

  double pulserH=0, pulserV=0, pulserDH=0, pulserDV=0; 

  TTree * friendly = new TTree("aux", "Auxdata (headers/gps)"); 
  friendly->Branch("header",&hdr); 
  friendly->Branch("pat",&patptr); 
  friendly->Branch("peakPulserCoherentH",&pulserH); 
  friendly->Branch("peakPulserCoherentV",&pulserV); 
  friendly->Branch("peakPulserDeconvolvedH",&pulserDH); 
  friendly->Branch("peakPulserDeconvolvedV",&pulserDV); 
  friendly->SetAutoSave(500); 

  for (size_t i = 0; i < strategies.size(); i++) 
  {
    TString tname; tname.Form("filter %s (on %s)",names[i].c_str(),label);
    trees[i] = new TTree(names[i].c_str(),tname.Data()); 
    trees[i]->Branch("summary",sum); 
    trees[i]->AddFriend(friendly); 
    trees[i]->SetAutoSave(500); 
  }


  time_t last_time = time(0); 
  bool auto_save = false;
  int ndone = 0; 

  for (int i = 0; i < d.N(); i++) 
  {
    if (time(0) - last_time > 1800)
    {
      auto_save = true; 
      last_time = time(0); 
    }

    d.getEntry(i); 

    hdr = d.header(); 
    if (isLDB) 
    {
      if (UCorrelator::isLDB(d.header()))
      {
         printf("----(LDB event %d ( idx=%d))-----\n",hdr->eventNumber,i); 
         FilteredAnitaEvent ev(d.useful(), &justAlfa, d.gps(), d.header()); 
         double phi,theta; 
         patptr->getThetaAndPhiWaveLDB(theta,phi);
         phi *= 180/TMath::Pi(); 
         theta *= 180/TMath::Pi(); 
         wc.combine(phi,theta,&ev,AnitaPol::kHorizontal); 
         pulserH = FFTtools::getPeakVal(wc.getCoherent()->hilbertEnvelope()); 
         pulserDH = FFTtools::getPeakVal(wc.getDeconvolved()->hilbertEnvelope()); 
         wc.combine(phi,theta,&ev,AnitaPol::kVertical); 
         pulserV = FFTtools::getPeakVal(wc.getCoherent()->hilbertEnvelope()); 
         pulserDV = FFTtools::getPeakVal(wc.getDeconvolved()->hilbertEnvelope()); 
         printf("  Pulser H: %g, Pulser V: %g\n", pulserH, pulserV); 

      }
      else continue; 
    }

    UsefulAdu5Pat pat(d.gps()); 


    if (isMC) 
    {
         printf("----(MC event %d ( idx=%d))-----\n",hdr->eventNumber,i); 
         FilteredAnitaEvent ev(d.useful(), &justAlfa, d.gps(), d.header()); 
         double phi,theta; 
         patptr->getThetaAndPhiWave(d.truth()->sourceLon, d.truth()->sourceLat, d.truth()->sourceAlt, theta,phi);
         phi *= 180/TMath::Pi(); 
         theta *= 180/TMath::Pi(); 
         wc.combine(phi,theta,&ev,AnitaPol::kHorizontal); 
         pulserH = FFTtools::getPeakVal(wc.getCoherent()->hilbertEnvelope()); 
         pulserDH = FFTtools::getPeakVal(wc.getDeconvolved()->hilbertEnvelope()); 
         wc.combine(phi,theta,&ev,AnitaPol::kVertical); 
         pulserV = FFTtools::getPeakVal(wc.getCoherent()->hilbertEnvelope()); 
         pulserDV = FFTtools::getPeakVal(wc.getDeconvolved()->hilbertEnvelope()); 
         printf("  Pulser H: %g, Pulser V: %g\n", pulserH, pulserV); 
    }

    if (isWAIS) 
    {
      if (UCorrelator::isWAISHPol(&pat, d.header()))
      {
         printf("----(WAIS event %d (idx=%d))-----\n",hdr->eventNumber,i); 
         FilteredAnitaEvent ev(d.useful(), &justAlfa, d.gps(), d.header()); 
         double phi,theta; 
         patptr->getThetaAndPhiWaveWaisDivide(theta,phi);
         phi *= 180/TMath::Pi(); 
         theta *= 180/TMath::Pi(); 
         wc.combine(phi,theta,&ev,AnitaPol::kHorizontal); 
         pulserH = FFTtools::getPeakVal(wc.getCoherent()->hilbertEnvelope()); 
         wc.combine(phi,theta,&ev,AnitaPol::kVertical); 
         pulserV = FFTtools::getPeakVal(wc.getCoherent()->hilbertEnvelope()); 
         printf("  Pulser H: %g, Pulser V: %g\n", pulserH, pulserV); 

      }
      else continue; 
    }

    if (isBG)
    {
      printf("----(event %d (idx=%d))-----\n",hdr->eventNumber,i); 
    }


    ofile.cd(); 
    patptr = &pat; 
    friendly->Fill(); 
    if (auto_save) friendly->AutoSave("SaveSelf"); 

    //preload all this stuff 
    d.useful(); 
    d.gps(); 
    d.header(); 

    for (size_t s = 0; s < strategies.size(); s++) 
    {
      printf (" %s...", names[s].c_str()); 
//  AnalysisWaveform::enableDebug(true); 
      FilteredAnitaEvent ev(d.useful(), strategies[s], d.gps(), d.header()); 
 // AnalysisWaveform::enableDebug(false); 
      analyzer->analyze(&ev, sum,d.truth()); 
      printf("[%g,%g]",sum->coherent[0][0].peakHilbert, sum->coherent[1][0].peakHilbert);
      ofile.cd(); 
      trees[s]->Fill(); 
      if (auto_save) trees[s]->AutoSave("SaveSelf"); 
    }
    printf("\n"); 

    if (max && ndone > max) break; 
    ndone++; 
  }


  ofile.cd(); 

  friendly->Write(); 
  for (size_t i = 0; i < trees.size(); i++) trees[i]->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
  ofile.Write(); 
}




