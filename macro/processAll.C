#include "FFTtools.h"
#include "TTree.h"
#include "AnitaDataset.h"
#include "TFile.h"
#include "TH2D.h"
#include "FilteredAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAnitaEvent.h"
#include "Adu5Pat.h"
#include "Analyzer.h"
#include "FilterStrategy.h"
#include "AnitaEventSummary.h"
#include "AnalysisConfig.h"
#include "BH13Filter.h"
#include "UCFilters.h"
#include "SystemResponse.h"
#include "TRandom3.h"

void processAll(int run)
{
	FFTtools::loadWisdom("wisdom.dat");
	AnitaDataset d(run);
	UCorrelator::AnalysisConfig cfg;

	cfg.nmaxima = 3;
  cfg.saturation_threshold = 1000;
  cfg.fill_channel_info = false;
	cfg.fill_blast_fraction = true;
  cfg.trace_to_continent = false;
	cfg.only_use_usable = true;
  cfg.zoomed_nant = 15;
  cfg.zoomed_ntheta = 30;
//  cfg.use_best_antenna_snr = true;

	cfg.response_option = UCorrelator::AnalysisConfig::ResponseA4;
	cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution;
	//cfg.fine_peak_finding_option = UCorrelator::AnalysisConfig::FinePeakFindingHistogram;
	
	UCorrelator::Analyzer analyzer(&cfg);
	double dtheta = 1.; double dphi = 2.; bool blockout = true;
	//analyzer.setTrackSun(dtheta, dphi, blockout);

	FilterStrategy forDeco;
	forDeco.addOperation(new UCorrelator::AntiBH13Filter());
	analyzer.setExtraFiltersDeconvolved(&forDeco);

	AnitaEventSummary* sum = new AnitaEventSummary;
	FilterStrategy strategy;
	UCorrelator::fillStrategyWithKey(&strategy, "sinsub_10_3_ad_2");
	strategy.addOperation(new UCorrelator::BH13Filter());

	TFile f(Form("/project2/avieregg/a4out/run%d.root", run), "RECREATE");
  TTree* tree = new TTree("sumTree", "sumTree");
	tree->Branch("summary", &sum);
/*
  TFile f(Form("/project2/avieregg/a4out/run%d.root", run), "UPDATE");
	TTree* tree = (TTree*) f.Get("sumTree");
	if(tree) tree->SetBranchAddress("summary", &sum);
  if(!tree)
  {
    tree = new TTree("sumTree", "sumTree");
    tree->Branch("summary", &sum);
  }
  tree->SetAutoSave(2000); 

  if(tree->GetEntries() == d.N())
  {
    fprintf(stderr, "run %d is finished, we're done here boys !!\n", run);
    return;
  }
*/
  d.getEntry(0);
  analyzer.getResponseManager()->checkTime(d.header()->payloadTime);

  int processed = 0;
	for(int i = tree->GetEntries(); i < d.N(); i++)
	//for(int i = 0; i < 100; i++)
	{
    if(i%1000 == 0) fprintf(stderr, "%d entries processed\n", i);
    //fprintf(stderr, "%d entries processed\n", i);
		d.getEntry(i);
		FilteredAnitaEvent fae(d.useful(), &strategy, d.gps(), d.header());

		analyzer.analyze(&fae, sum);
		f.cd();
		tree->Fill();
    processed++;
    //if(processed == 30000) break;
	}
  fprintf(stderr, "run %d is finished, we're done here boys !!\n", run);
	f.cd();
	tree->Write();
	f.Close();
	FFTtools::saveWisdom("wisdom.dat");
}
