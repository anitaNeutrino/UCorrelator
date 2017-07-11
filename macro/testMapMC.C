#include "ProbabilityMap.h" 
#include "UCUtil.h"
#include "UCFilters.h" 
#include "PointingResolutionModel.h"
#include "Analyzer.h" 
#include "AnalysisConfig.h" 
#include "AnitaDataset.h" 
#include "SystemResponse.h" 


UCorrelator::ProbabilityMap* testMapMC(int run =223, int max = 1000) 
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

  StereographicGrid g(4096,4096); 
  UCorrelator::ConstantPointingResolutionModel m(0.2,0.3);
  UCorrelator::ProbabilityMap *map = new UCorrelator::ProbabilityMap(&g,&m); 
  int ndone = 0; 

  int eventNumber = 0;
  for (int i = 0; i< d.N(); i++)
  {
    d.getEntry(i); 
    if (d.header()->eventNumber == eventNumber) continue; 
    printf("----(%d, %d)-----\n",i, d.header()->eventNumber); 
    eventNumber = d.header()->eventNumber; 

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

  
  map->segmentationScheme()->Draw("colz", map->getProbabilities()); 

  TFile f("g.root"); 
  TGraph * gg = (TGraph*) f.Get("Graph"); 
  gg->SetLineColor(1); 
  gg->Draw("lsame"); 

  return map; 
}
