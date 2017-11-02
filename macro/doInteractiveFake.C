#include "FFTtools.h"

UCorrelator::Analyzer *doInteractiveFake(int run = 342, int event = 0, bool simulated = false )
{

  FFTtools::loadWisdom("wisdom.dat"); 
  FilterStrategy *strategy= new FilterStrategy; 
 

  AnitaDataset d(run,false, WaveCalType::kDefault, simulated ? AnitaDataset::ANITA_MC_DATA : AnitaDataset::ANITA_ROOT_DATA, AnitaDataset::kNoBlinding );

  event > 0 ? d.getEvent(event) : d.getEntry(-event); 



  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 1; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
  cfg.max_peak_trigger_angle =60; 
  cfg.use_forced_trigger_rms = false; 
//  cfg.correlator_theta_lowest=90; 
//  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 
//  cfg.response_option = UCorrelator::AnalysisConfig::ResponseHarmSignalOnly; 
  //cfg.combine_unfiltered = false; 

  UCorrelator::Analyzer * analyzer = new UCorrelator::Analyzer(&cfg,true); 

//  strategy->addOperation(new UCorrelator::DeconvolveFilter(analyzer->getResponseManager(), cfg.deconvolution_method)); 
//  strategy->addOperation(new ALFAFilter); 
//  strategy->addOperation(ssf); 

//  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("adsinsub_3_10_3"), d.gps(), d.header()); 
  AnitaEventFaker faker("IndividualBRotter"); 
  faker.makePureNoiseEvent(0.1,d.useful()); 
  faker.addSignal(d.useful(),10,175,10, std::complex<double>(1/sqrt(2),0), std::complex<double>(-1/sqrt(2),0)); 
  FilteredAnitaEvent* ev = new FilteredAnitaEvent(d.useful(),UCorrelator::getStrategyWithKey("sinsub_10_3_ad_2"), d.gps(), d.header()); 
  printf("auto fae = (FilteredAnitaEvent *) %p;\n",ev); 

//  ev->plotSummary(0,0); 


  AnitaEventSummary sum; 
  analyzer->analyze(ev,&sum,d.truth()); 
  analyzer->drawSummary(); 
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
  TCanvas * filter =  new TCanvas("filter","Filter"); 
  filter->Divide(2,1); 
  filter->cd(1); 
  mp->getCurrentFilterPower(AnitaPol::kHorizontal,0)->Draw(); 
  filter->cd(2); 
  mp->getCurrentFilterTimeDomain(AnitaPol::kHorizontal,0)->Draw(); 
  */

   

//  butter->getFilter(AnitaPol::kHorizontal,0)->drawResponse(0,101,10); 

  FFTtools::saveWisdom("wisdom.dat"); 
  return analyzer; 
}
