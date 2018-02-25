
#include "FFTtools.h"

UCorrelator::ProbabilityMap::Params * map_params()
{
  // pixel x, pixel y , max meter x, max meter y
  StereographicGrid * g= new StereographicGrid(1000,1000,2000000,2000000); 

  TF1 * f_dtheta = new TF1("ftheta", "[0] / x^[1]", 1, 50);
  TF1 * f_dphi = new TF1("fphi", "[0] / x^[1]", 1, 50);
  //anita4 fit from wais.
  f_dtheta->SetParameter(0, 5.431); 
  f_dtheta->SetParameter(1, 1.155); 
  f_dphi->SetParameter(0, 28.87); 
  f_dphi->SetParameter(1, 1.398);
  //anita3
  // f_dtheta->SetParameter(0, 0.3936); 
  // f_dtheta->SetParameter(1, 0.2102); 
  // f_dphi->SetParameter(0, 1.065); 
  // f_dphi->SetParameter(1, 0.2479); 
  UCorrelator::PointingResolutionParSNRModel * snrResolutionModel = new UCorrelator::PointingResolutionParSNRModel (*f_dtheta, *f_dphi, true,true);
  UCorrelator::PointingResolutionModelPlusHeadingError * resolutionModel = new UCorrelator::PointingResolutionModelPlusHeadingError(20, snrResolutionModel); 

  Refraction::SphRay * ref = new Refraction::SphRay; 

  UCorrelator::ProbabilityMap::Params *p = new UCorrelator::ProbabilityMap::Params; 
  p->refract = ref; 
  // p->refract = 0; 
  p->seg = g; 
  p->point = resolutionModel; 
  p->collision_detection = true; 
  p->verbosity = 0; // verbosity level for output info.
  p->maximum_distance = 2.5;
  p->min_p_on_continent = 0;
  return p; 

}

// UCorrelator::Analyzer *doInteractive(int run = 130, int event = 22896140, bool decimated = false, bool simulated = false )
// UCorrelator::Analyzer *doInteractive(int run = 140, int event = 25639095, bool decimated = false, bool simulated = false )
// UCorrelator::Analyzer *doInteractive(int run = 153, int event = 30003847, bool decimated = false, bool simulated = false )
// UCorrelator::Analyzer *doInteractive(int run = 102, int event = 204382, bool decimated = false, bool simulated = true )
// UCorrelator::Analyzer *doInteractive(int run = 102, int event = 204407, bool decimated = false, bool simulated = true )
// UCorrelator::Analyzer *doInteractive(int run = 160, int event = 32096871, bool decimated = false, bool simulated = false )
// UCorrelator::Analyzer *doInteractive(int run = 308, int event = 80224937, bool decimated = false, bool simulated = false )
// UCorrelator::Analyzer *doInteractive(int run = 349, int event = 92484106, bool decimated = false, bool simulated = false )
UCorrelator::Analyzer *doInteractive(int run = 3, int event = 7579449, bool decimated = false, bool simulated = true )
// UCorrelator::Analyzer *doInteractive(int run = 7, int event = 14305756, bool decimated = false, bool simulated = true )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy *strategy= new FilterStrategy; 
 
//  TF1 fn("foo"," (x < 0.2) * exp((x-0.2)/0.01)  + (x > 0.2 && x < 1.2) * (1-0.05*x) + (x > 1.2) * exp((1.2-x)/0.02)", 0,2); 

//  UCorrelator::SpectrumAverage *avg = new UCorrelator::SpectrumAverage(run,60); 
//  avg->computePeakiness(); 

 // UCorrelator::CombinedSineSubtractFilter * ssf = new UCorrelator::CombinedSineSubtractFilter(0.05,10);
 // ssf->setInteractive(true); 
//  printf("UCorrelator::CombinedSineSubtractFilter* ssf = (UCorrelator::CombinedSineSubtractFilter*) %p\n", ssf); 
//  strategy->addOperation(ssf);
//  strategy->addOperation(new SimplePassBandFilter(0.2,1.2)); 

// UCorrelator::SineSubtractFilter * ssf = new UCorrelator::SineSubtractFilter(0.10,3);
//  ssf->makeAdaptive(avg); 
//  UCorrelator::AdaptiveMinimumPhaseFilter * mp = new UCorrelator::AdaptiveMinimumPhaseFilter(avg,-2,5); 
//  printf("UCorrelator::AdaptiveMinimumPhaseFilter * mp = (UCorrelator::AdaptiveMinimumPhaseFilter *) %p\n",mp); 
//  strategy->addOperation(mp); 
//
//UCorrelator::AdaptiveButterworthFilter * butter = new UCorrelator::AdaptiveButterworthFilter(&avg); 
//printf("UCorrelator::AdaptiveButterworthFilter * butter = (UCorrelator::AdaptiveButterworthFilter *) %p\n",butter); 
//strategy->addOperation(butter); 


  AnitaDataset d(run,decimated, WaveCalType::kDefault, simulated ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA, AnitaDataset::kNoBlinding );

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 



 UCorrelator::AnalysisConfig cfg;
    cfg.nmaxima = 3;
    cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;
    cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;
  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 

  FilterStrategy* forDeco = new FilterStrategy;
  forDeco->addOperation(new UCorrelator::AntiBH13Filter());
  analyzer->setExtraFiltersDeconvolved(forDeco);
  analyzer->setDisallowedAntennas(0, (1ul<<45));  // Vpol ant45 is bad! So disable it.

  UCorrelator::fillStrategyWithKey(strategy, "sinsub_10_3_ad_2");
  strategy->addOperation(new UCorrelator::BH13Filter()); 

//  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("adsinsub_3_10_3"), d.gps(), d.header()); 
  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),strategy, d.gps(), d.header()); 
  printf("auto fae = (FilteredAnitaEvent *) %p;\n",ev); 

  // ev->plotSummary(0,0); 


  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum,d.truth()); 
  analyzer->drawSummary(); 
  std::cout<< "eventNumber "<<sum.eventNumber<<std::endl;
  /*
  TCanvas * c2 = new TCanvas; 
  c2->Divide(2,1); 
  const UCorrelator::AbstractResponse * response = analyzer->getResponseManager()->response(0,0); 
  AnalysisWaveform * imp =response->impulseResponse(0.1, 981); 
  c2->cd(1); 
  imp->drawEven(); 
  AnalysisWaveform * deconv = response->deconvolve(imp,cfg.deconvolution_method); 
  c2->cd(2); 
  deconv->drawEven(); 
  */

  /*
  TCanvas * "sinsub_10_3_ad_2" =  new TCanvas(""sinsub_10_3_ad_2"",""sinsub_10_3_ad_2""); 
  "sinsub_10_3_ad_2"->Divide(2,1); 
  "sinsub_10_3_ad_2"->cd(1); 
  mp->getCurrentFilterPower(AnitaPol::kHorizontal,0)->Draw(); 
  "sinsub_10_3_ad_2"->cd(2); 
  mp->getCurrentFilterTimeDo  main(AnitaPol::kHorizontal,0)->Draw(); 
  */

  TCanvas * psCanvas = new TCanvas();
  psCanvas->cd();
  double p_ground;
  UCorrelator::ProbabilityMap::Params * p = map_params(); 
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(p); 
  map->add(p_ground, &sum, d.gps(), AnitaPol::AnitaPol_t(sum.mostImpulsivePolAsInt()), sum.mostImpulsiveInd(), 1);
  map->segmentationScheme()->Draw("colz",map->getProbSums(false));




//  butter->getFilter(AnitaPol::kHorizontal,0)->drawResponse(0,101,10); 

  FFTtools::saveWisdom("wisdom.dat"); 
  return analyzer; 
}
