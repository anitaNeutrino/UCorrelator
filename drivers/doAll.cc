#include "FFTtools.h"
#include "AnitaConventions.h"
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
#include "Hical2.h"


void doAll( int run = 352, int max = 0, int start = 0, const char * filter = "sinsub_10_3_ad_2" )
//void doWais2( int run = 352, int max = 0, int start = 0, const char * filter = "" )
{

  FFTtools::loadWisdom("wisdom.dat");

//  /*AnalysisWaveform::InterpolationType*/ AnalysisWaveform::defaultInterpolationType = AnalysisWaveform::REGULARIZED_SPARSE_YEN;

  AnitaDataset d(run,false,WaveCalType::kDefault, AnitaDataset::ANITA_ROOT_DATA, AnitaDataset::kNoBlinding );
  UCorrelator::AnalysisConfig cfg;
    cfg.nmaxima = 3;
    cfg.response_option = UCorrelator::AnalysisConfig::ResponseTUFF;
    cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;
  UCorrelator::Analyzer analyzer(&cfg);

  TString outname;
  if (max && start) outname.Form("a%dall/%d_max_%d_start_%d_%s.root",AnitaVersion::get(),run,max,start,filter);
  else if (max) outname.Form("a%dall/%d_max_%d_%s.root",AnitaVersion::get(),run,max,filter);
  else if (start) outname.Form("a%dall/%d_start_%d_%s.root",AnitaVersion::get(),run,start,filter);
  else outname.Form("a%dall/%d_%s.root",AnitaVersion::get(),run, filter);

  TFile ofile(outname, "RECREATE");
  TTree tree(TString::Format("anita%d",AnitaVersion::get()),TString::Format("anita%d", AnitaVersion::get())); 
  tree.SetAutoFlush(1000); 
  AnitaEventSummary * sum = new AnitaEventSummary();


  double dtheta = 1.; double dphi = 2.; bool blockout = true;
  analyzer.setTrackSun(dtheta, dphi, blockout);

  FilterStrategy* forDeco = new FilterStrategy;
  forDeco->addOperation(new UCorrelator::AntiBH13Filter());
  analyzer.setExtraFiltersDeconvolved(forDeco);
  analyzer.setDisallowedAntennas(0, (1ul<<45));  // Vpol ant45 is bad! So disable it.

  FilterStrategy strategy (&ofile);
  UCorrelator::fillStrategyWithKey(&strategy, filter);
  strategy.addOperation(new UCorrelator::BH13Filter());


  RawAnitaHeader *hdr = 0 ;
  UsefulAdu5Pat *patptr = 0;
  // double isHC = 0;
  tree.Branch("summary",&sum);
  tree.Branch("header",&hdr);
  tree.Branch("pat",&patptr);
  // tree.Branch("isHC",&isHC);

  int ndone = 0;
  for (int i =start ; i < d.N(); i++)
  {
      d.getEntry(i);
      printf("----(%d)-----\n",i);

      UsefulAdu5Pat pat(d.gps());
      // 1% data, except the wais events
      if (TString::Hash(&d.header()->eventNumber, sizeof(d.header()->eventNumber))%10 == 0 && !UCorrelator::isWAISHPol(&pat, d.header()) && !UCorrelator::isWAISVPol(&pat, d.header()))
      {
        const time_t ctt = time(0);
        printf("Processing event %d (%d) \t|%s", d.header()->eventNumber,ndone,asctime(localtime(&ctt)));        
        FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header());

        analyzer.analyze(&ev, sum);
        ofile.cd();
        hdr = d.header();
        patptr = &pat;
        // isHC = Hical2::isHical(sum);
        tree.Fill();
        ndone++;
      }
      if (max && ndone >= max) break;
  }

  ofile.cd();
  tree.Write();

  FFTtools::saveWisdom("wisdom.dat");
  std::cout << "end of script"<<std::endl;
}

int main (int nargs, char ** args)
{

  int run = nargs < 2 ? 352 : atoi(args[1]);
  int max = nargs < 3 ? 0 : atoi(args[2]);
  int start = nargs < 4 ? 0 : atoi(args[3]);
  const char * filter = nargs < 5 ? 0 :args[4];

  if (filter)
    doAll(run,max,start,filter);
  else
    doAll(run,max,start);
}