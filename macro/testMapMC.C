#include "ProbabilityMap.h" 
#include "UCUtil.h"
#include "UCFilters.h" 
#include "PointingResolutionModel.h"
#include "Analyzer.h" 
#include "AnalysisConfig.h" 
#include "AnitaDataset.h" 
#include "SystemResponse.h" 


UCorrelator::ProbabilityMap* testMapMC(int run =223, int max = 10) 
{

  FFTtools::loadWisdom("wisdom.dat"); 
  AnitaDataset d(run,false,WaveCalType::kDefault, AnitaDataset::ANITA_MC_DATA); 


  UCorrelator::AnalysisConfig cfg; 
  cfg.nmaxima = 2; 
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseIndividualBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution; 

  UCorrelator::Analyzer  analyzer (&cfg); 
  AnitaEventSummary sum; 
  FilterStrategy strategy; 
  UCorrelator::fillStrategyWithKey(&strategy,"sinsub_10_3_ad_2"); 

  StereographicGrid g(1024,1024); 
  UCorrelator::ConstantPointingResolutionModel m(0.2,0.3);
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&g,&m); 
  int ndone = 0; 

  for (int i = 0; i< d.N(); i++)
  {
    d.getEntry(i); 
    printf("----(%d)-----\n",i); 
    UsefulAdu5Pat pat(d.gps()); 
    FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header()); 
    analyzer.analyze(&ev, &sum,d.truth()); 


    int ipol = sum.mc.wf[0].peakHilbert > sum.mc.wf[1].peakHilbert ? 0 : 1; 

    for (int j = 0; j < cfg.nmaxima;j++)
    {

      if ( fabs(FFTtools::wrap(sum.mc.phi-sum.peak[ipol][j].phi,360,0)) < 3 
          && fabs(FFTtools::wrap(sum.mc.theta-sum.peak[ipol][j].theta,360,0)) < 2)
      {

        map->add(&sum,d.gps(), AnitaPol::AnitaPol_t (ipol), j); 
        ndone++;
      }
    }

    if (max && ndone >= max) break; 
  }


  FFTtools::saveWisdom("wisdom.dat"); 

  AntarcticaBackground * bg = new AntarcticaBackground(); 
  bg->Draw("colz"); 
  
  map->segmentationScheme()->Draw("colzsame", map->getProbabilities()); 
  return map; 
}
