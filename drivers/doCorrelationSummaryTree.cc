//////////////////////////////////////////////////////
//
//  Program to make CorrelationSummaryTrees
//    with different filter strategies
//
//////////////////////////////////////////////////////

#include "FFTtools.h"
#include "Analyzer.h"
#include "FilteredAnitaEvent.h"
#include "PrettyAnitaEvent.h"
#include "AnitaEventCalibrator.h"
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
#include "CorrelationSummaryAnita3.h"

void doCorrelationSummaryTree( int run = 352, int max = 0, int start = 0, const char * filter = "sinsub_5_3_ad_2", const char * outputDir = "correlationSummary")
{
  
  FFTtools::loadWisdom("wisdom.dat"); 

  AnitaVersion::set(3);

  AnitaGeomTool *fGeomTool = AnitaGeomTool::Instance();
  fGeomTool->usePhotogrammetryNumbers(1);
  AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  for(Int_t surf=0; surf<NUM_SURF; surf++){
    for(Int_t chan=0; chan<NUM_CHAN; chan++){
      cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0; ///< From phase center to AMPAs (hopefully)
    }
  }

  Double_t deltaT= 1. / (2.6*40.);
  
  AnitaPol::AnitaPol_t pol;

  Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
  Int_t cutTimeNs = 1200;//60e3;
  char cpol[100];
  bool isLDB = false;

  Int_t delayGenerator = 199996850; // delay generator was not exactly 200 ms
  Int_t delay=0;
  Int_t constdelay = 500;
  Int_t deltaTriggerTimeNs;

  if (run>300){ // WAIS
    pol = AnitaPol::kHorizontal;
  } else if (run<150){ // LDB VPOL 
    pol = AnitaPol::kVertical;
    isLDB = true;
    delay =  25000000; // V-POL pulse at 25 ms runs 145-149
  } else if (run<154){ // LDB HPOL
    pol = AnitaPol::kHorizontal;
    isLDB = true;
    delay =  50000000; // H-POL pulse at 50 ms
  } else if (run<172){ // LDB VPOL
    pol = AnitaPol::kVertical;
    isLDB = true;
    delay =  50000000; // V-POL pulse at 50 ms runs 154-171
  } else {
    std::cout << "Unknown run" << std::endl;
    return;
  }

  if (pol == AnitaPol::kVertical) sprintf(cpol, "VPOL");
  else if (pol == AnitaPol::kHorizontal) sprintf(cpol, "HPOL");
  
  
  AnitaDataset d(run); 
  UCorrelator::AnalysisConfig cfg; 
  cfg.start_pol = AnitaPol::kHorizontal; 
  cfg.end_pol = AnitaPol::kHorizontal; 
  
  cfg.response_option = UCorrelator::AnalysisConfig::ResponseSingleBRotter; 
  cfg.deconvolution_method = new AnitaResponse::AllPassDeconvolution(); 


  UCorrelator::Analyzer analyzer(&cfg); 


  
  TString outname; 
  outname.Form("%s/CorrelationSummaryTree_%d_%s.root", outputDir, run, filter); 

  TFile ofile(outname, "RECREATE"); 

  FilterStrategy strategy (&ofile); 
  UCorrelator::fillStrategyWithKey(&strategy, filter); 

  
  CorrelationSummaryAnita3 *theCor=0;
  TTree *corTree = new TTree("corTree","Tree of Correlation Summaries");
  corTree->Branch("cor","CorrelationSummaryAnita3",&theCor);

  FilteredAnitaEvent *ev = 0;
  Double_t thetaWave, phiWave;
  Int_t ant;
  int ndone=0;
  
  for (int i =start ; i < d.N(); i++)
  {

    d.getEntry(i); 
    //    printf("----(%d)-----\n",i); 

    UsefulAdu5Pat pat(d.gps()); 

    if (UCorrelator::isWAISHPol(&pat, d.header()))
    {
      printf("Processing event %d (%d)\n",d.header()->eventNumber,ndone); 
      ev = new FilteredAnitaEvent(d.useful(), &strategy, d.gps(), d.header());      

      pat.getThetaAndPhiWaveWaisDivide(thetaWave,phiWave);
      
      ant=fGeomTool->getTopAntNearestPhiWave(phiWave, pol);

      PrettyAnitaEvent dumb(d.calibrated());
      theCor = dumb.createCorrelationSummaryAnita3(ant,pol,deltaT);

      for(int corInd=0;corInd<NUM_CORRELATIONS_ANITA3;corInd++) {

	const TGraphAligned *g1 = ev->getFilteredGraph(theCor->firstAnt[corInd],  pol)->even();
	const TGraphAligned *g2 = ev->getFilteredGraph(theCor->secondAnt[corInd], pol)->even();

	TGraph *grCor;

	if(deltaT==0) {
	  grCor = FFTtools::getCorrelationGraph(g1, g2);
	}
	else {
	  grCor = FFTtools::getInterpolatedCorrelationGraph(g1, g2, deltaT);
	}

	double *theTimes = grCor->GetX();
	double *theValues = grCor->GetY();
      
	int numPoints=grCor->GetN();

	double rmsVal=TMath::RMS(numPoints,theValues);
	int maxIndex=TMath::LocMax(numPoints,theValues);

	theCor->rmsCorVals[corInd]=rmsVal;
	theCor->maxCorVals[corInd]=theValues[maxIndex];
	theCor->maxCorTimes[corInd]=theTimes[maxIndex];

	// std::cout << theCor->firstAnt[corInd] << "\t" << theCor->secondAnt[corInd]
	// 	  << "\t" << theCor->maxCorTimes[corInd] 
	// 	  << "\t" << theCor->maxCorVals[corInd] << "\t" 
	// 	  << "\t" << std::endl;

	theCor->secondCorVals[corInd][0]=theCor->maxCorVals[corInd];
	theCor->secondCorTimes[corInd][0]=theCor->maxCorTimes[corInd];
	theCor->secondCorVals[corInd][1]=theCor->maxCorVals[corInd];
	theCor->secondCorTimes[corInd][1]=theCor->maxCorTimes[corInd];
	for(int i=maxIndex-1;i>=1;i--) {
	  if(i<1) break;	 
	  if(theValues[i]>=theValues[i-1] && theValues[i]>=theValues[i+1]) {
	    theCor->secondCorVals[corInd][0]=theValues[i];
	    theCor->secondCorTimes[corInd][0]=theTimes[i];
	    break;
	  }	  
	}
	for(int i=maxIndex+1;i<grCor->GetN();i++) {
	  if(i>=grCor->GetN()-1) break;	 
	  if(theValues[i]>=theValues[i-1] && theValues[i]>=theValues[i+1]) {
	    theCor->secondCorVals[corInd][1]=theValues[i];
	    theCor->secondCorTimes[corInd][1]=theTimes[i];
	    break;
	  }	  
	}
	
	delete grCor;
      }
      
      
      corTree->Fill();     
      delete theCor;
      delete ev;     
      
      ofile.cd(); 

      ndone++;

      
    }

    if (max && ndone >= max) break; 

  }

  ofile.cd(); 
  corTree->Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
}

int main (int nargs, char ** args)
{
   
  int run = nargs < 2 ? 352 : atoi(args[1]);
  int max = nargs < 3 ? 0 : atoi(args[2]); 
  int start = nargs < 4 ? 0 : atoi(args[3]); 
  const char * filter = nargs < 5 ? 0 :args[4]; 
  const char * outDir = nargs < 6 ? 0 :args[5]; 

  if (filter) {
    if (outDir)
      doCorrelationSummaryTree(run, max, start, filter, outDir);
    else
      doCorrelationSummaryTree(run, max, start, filter);
  }
  else
    doCorrelationSummaryTree(run, max, start);


}
