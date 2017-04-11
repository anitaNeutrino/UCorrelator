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

void addStrategy( FilterStrategy * s, const char * name) { strategies.push_back(s); names.push_back(name); } 
UCorrelator::Analyzer * analyzer; 
UCorrelator::AnalysisConfig cfg; 


/** add filter strategies here ! */ 

void setupFilters(TFile* out, int run) 
{
  (void) out; 
  justAlfa.addOperation(new ALFAFilter); 

  /** Sine subtraction */ 
  FilterStrategy * sinsub_05_0= new FilterStrategy; 
  UCorrelator::SineSubtractFilter * ssf = new UCorrelator::SineSubtractFilter(0.05, 0); 
  sinsub_05_0->addOperation(ssf); 
  sinsub_05_0->addOperation(new ALFAFilter); 
  addStrategy(sinsub_05_0, "sinsub_05_0"); 

  FilterStrategy * sinsub_03_0= new FilterStrategy; 
  ssf = new UCorrelator::SineSubtractFilter(0.03, 0); 
  sinsub_03_0->addOperation(ssf); 
  sinsub_03_0->addOperation(new ALFAFilter); 
  addStrategy(sinsub_03_0, "sinsub_03_0"); 


  FilterStrategy * sinsub_10_0= new FilterStrategy; 
  sinsub_10_0->addOperation(new UCorrelator::SineSubtractFilter(0.1, 0)); 
  sinsub_10_0->addOperation(new ALFAFilter); 
  addStrategy(sinsub_10_0, "sinsub_10_0"); 

  FilterStrategy * sinsub_05_3= new FilterStrategy; 
  sinsub_05_3->addOperation(new UCorrelator::SineSubtractFilter(0.05, 3)); 
  sinsub_05_3->addOperation(new ALFAFilter); 
  addStrategy(sinsub_05_3, "sinsub_05_3"); 


  UCorrelator::SpectrumAverage *  avg  = new UCorrelator::SpectrumAverage(run,60,"specavg"); //TODO use default dir for this 
  avg->computePeakiness(); 

  /** Peakiness Detecting Sine Subtraction */ 

  FilterStrategy * adsinsub_05_0= new FilterStrategy; 
  UCorrelator::SineSubtractFilter * adssf = new UCorrelator::SineSubtractFilter(0.05, 0); 
  adssf->makeAdaptive(avg); 
  adsinsub_05_0->addOperation(adssf); 
//  adssf->setVerbose(1); 
  adsinsub_05_0->addOperation(new ALFAFilter); 
  addStrategy(adsinsub_05_0, "adsinsub_05_0"); 

  FilterStrategy * adsinsub_10_0= new FilterStrategy; 
  adssf = new UCorrelator::SineSubtractFilter(0.10, 0); 
  adssf->makeAdaptive(avg); 
//  adssf->setVerbose(1); 
  adsinsub_10_0->addOperation(adssf); 
  adsinsub_10_0->addOperation(new ALFAFilter); 
  addStrategy(adsinsub_10_0, "adsinsub_10_0"); 

  FilterStrategy * adsinsub_10_3= new FilterStrategy; 
  adssf = new UCorrelator::SineSubtractFilter(0.10, 3); 
  adssf->makeAdaptive(avg); 
//  adssf->setVerbose(1); 
  adsinsub_10_3->addOperation(adssf); 
  adsinsub_10_3->addOperation(new ALFAFilter); 
  addStrategy(adsinsub_10_3, "adsinsub_10_3"); 



  /** deconv + sine subtraction */ 
  FilterStrategy * decon_ss_5_0 = new FilterStrategy; 
  decon_ss_5_0->addOperation(new UCorrelator::DeconvolveFilter(analyzer->getResponseManager(), cfg.deconvolution_method)); 
  decon_ss_5_0->addOperation(new ALFAFilter); 
  decon_ss_5_0->addOperation(new UCorrelator::SineSubtractFilter(0.05,0)); 
  addStrategy(decon_ss_5_0,"decon_ss_5_0"); 

  FilterStrategy * decon_ss_3_0 = new FilterStrategy; 
  decon_ss_3_0->addOperation(new UCorrelator::DeconvolveFilter(analyzer->getResponseManager(), cfg.deconvolution_method)); 
  decon_ss_3_0->addOperation(new ALFAFilter); 
  decon_ss_3_0->addOperation(new UCorrelator::SineSubtractFilter(0.03,0)); 
  addStrategy(decon_ss_3_0,"decon_ss_3_0"); 

  FilterStrategy * decon_adss_5_0 = new FilterStrategy; 
  decon_adss_5_0->addOperation(new UCorrelator::DeconvolveFilter(analyzer->getResponseManager(), cfg.deconvolution_method)); 
  decon_adss_5_0->addOperation(new ALFAFilter); 
  adssf = new UCorrelator::SineSubtractFilter(0.05,0); 
  adssf->makeAdaptive(avg); 
  decon_adss_5_0->addOperation(adssf); 
 addStrategy(decon_adss_5_0,"decon_adss_5_0"); 


  /** Adaptive Butterworth Filter */ 
  FilterStrategy * butter_2 = new  FilterStrategy; 
  butter_2->addOperation(new UCorrelator::AdaptiveButterworthFilter(avg,2)); 
  butter_2->addOperation(new ALFAFilter); 
  addStrategy(butter_2, "butter_2"); 

  FilterStrategy * butter_15 = new  FilterStrategy; 
  butter_15->addOperation(new UCorrelator::AdaptiveButterworthFilter(avg,1.5)); 
  butter_15->addOperation(new ALFAFilter); 
  addStrategy(butter_15, "butter_15"); 



  /** Adaptive Minimum Phase exponent = -2**/ 
  FilterStrategy * minphase_2 = new  FilterStrategy; 
  minphase_2->addOperation(new UCorrelator::AdaptiveMinimumPhaseFilter(avg,-2)); 
  minphase_2->addOperation(new ALFAFilter); 
  addStrategy(minphase_2, "minphase_2"); 

  FilterStrategy * minphase_1 = new  FilterStrategy; 
  minphase_1->addOperation(new UCorrelator::AdaptiveMinimumPhaseFilter(avg,-1)); 
  minphase_1->addOperation(new ALFAFilter); 
  addStrategy(minphase_1, "minphase_1"); 

  /**Geometric Filter */ 

  FilterStrategy * geom = new FilterStrategy; 
  geom->addOperation(new ALFAFilter); 

  std::vector<std::vector<TGraphAligned* > > noise(48, std::vector<TGraphAligned*>(2)); 
  for (int ant = 0; ant < 48; ant++)
  {
    for (int pol = 0; pol < 2; pol++) 
    {
      TH1 * avg = UCorrelator::SpectrumAverage::defaultThermal()->getSpectrumPercentile(AnitaPol::AnitaPol_t(pol), ant,0.5,true); 
      noise[ant][pol] = new TGraphAligned(avg->GetNbinsX()); 
      for (int i = 0; i < avg->GetNbinsX(); i++) 
      {
        noise[ant][pol]->GetX()[i] = avg->GetBinLowEdge(i+1); 
        noise[ant][pol]->GetY()[i] = avg->GetBinContent(i+1); 
      }
    }
  }

  FilterOperation * geomfilter = new GeometricFilter(noise); 
  geom->addOperation(geomfilter); 

  addStrategy(geom,"geom"); 

}




int main(int nargs, char ** args) 
{

  AnitaVersion::set(3); 
  FFTtools::loadWisdom("wisdom.dat"); 
  int run = atoi(args[1]); 

  bool isWAIS = run >=min_wais && run <= max_wais; 
  bool isLDB = run >=min_ldb && run <= max_ldb; 
  bool isBG  = !isWAIS && !isLDB; 

  AnitaDataset d(run, isBG); //only use decimated if background 

  int max = nargs > 2 ? atoi(args[2]) : 0; 

  TString outname; 
  const char * label = isWAIS ? "wais" : isLDB ? "ldb" : "bg"; 

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



  setupFilters(&ofile, run); 
 
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
      analyzer->analyze(&ev, sum); 
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




