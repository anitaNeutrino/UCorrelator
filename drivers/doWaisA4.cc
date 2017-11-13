//#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
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
#include "BH13Filter.h"


void doWaisA4( int run = 352, int max = 0, int start = 0, const char * filter = "sinsub_10_3_ad_2" )
//void doWais2( int run = 352, int max = 0, int start = 0, const char * filter = "" )
{

  FFTtools::loadWisdom("wisdom.dat");

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN;

  AnitaDataset d(run,false,WaveCalType::kDefault, AnitaDataset::ANITA_ROOT_DATA, AnitaDataset::kNoBlinding );
  UCorrelator::AnalysisConfig cfg;
    //cfg.nmaxima = 3;
    cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;
    //cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;
  UCorrelator::Analyzer analyzer(&cfg);

  TString outname;
  if (max && start) outname.Form("wais/%d_max_%d_start_%d_%s.root",run,max,start,filter);
  else if (max) outname.Form("wais/%d_max_%d_%s.root",run,max,filter);
  else if (start) outname.Form("wais/%d_start_%d_%s.root",run,start,filter);
  else outname.Form("wais/%d_%s.root",run, filter);

  TFile ofile(outname, "RECREATE");
  TTree * tree = new TTree("wais","WAIS Hpol");
  AnitaEventSummary * sum = new AnitaEventSummary;


    double dtheta = 1.; double dphi = 2.; bool blockout = true;
    analyzer.setTrackSun(dtheta, dphi, blockout);

    FilterStrategy* forDeco = new FilterStrategy;
    forDeco->addOperation(new UCorrelator::AntiBH13Filter());
    analyzer.setExtraFiltersDeconvolved(forDeco);

  FilterStrategy strategy (&ofile);
  UCorrelator::fillStrategyWithKey(&strategy, filter);
  strategy.addOperation(new UCorrelator::BH13Filter());


  RawAnitaHeader *hdr = 0 ;
  UsefulAdu5Pat *patptr = 0;
  tree->Branch("summary",&sum);
  tree->Branch("header",&hdr);
  tree->Branch("pat",&patptr);

  int ndone = 0;
  for (int i =start ; i < d.N(); i++)
  {
    try{
      d.getEntry(i);
      printf("----(%d)-----\n",i);

      UsefulAdu5Pat pat(d.gps());

      if (UCorrelator::isWAISHPol(&pat, d.header()) || UCorrelator::isWAISVPol(&pat, d.header()))
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

    }catch(const char* msg){
      std::cout<<"an error catched for this event!!!!"<< msg <<std::endl;
    }

    

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

  if (filter)
    doWais2(run,max,start,filter);
  else
    doWais2(run,max,start);


}
