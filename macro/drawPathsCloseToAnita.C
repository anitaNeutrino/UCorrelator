#include "TFile.h"
#include "TTree.h"
#include "AnitaEventSummary.h" 
#include "FFTtools.h"
#include "BaseList.h"
// #include "UsefulAdu5Pat.h"
// .x drawPathsCloseToAnita.C
//copy the top headers in the root console and run .x drawPathsCloseToAnita.C
void drawPathsCloseToAnita(const TString fileName = "/Users/sylarcp/anitaNeutrino/anitaBuildTool/components/UCorrelator/macro/sparsedAllRuns/sparsedAllRuns.root"){
  AnitaEventSummary * sum = new AnitaEventSummary; 
  Adu5Pat * pat = new Adu5Pat; 
  TFile * sumfile = new TFile(fileName); 
  TTree * sumtree = (TTree*) sumfile->Get("anita4"); 
  sumtree->SetBranchAddress("summary",&sum); 
  sumtree->SetBranchAddress("pat",&pat); 
  int N=sumtree->GetEntries();
  const double maxDistKm = 800;
  // const double maxDistKm = 8000000;
  BaseList::makeBaseList();
  std::vector<TGraphAntarctica*> grs(BaseList::getNumAbstractBases(), NULL);
  TGraphAntarctica* grFlightPath = new TGraphAntarctica();
  grFlightPath->SetName("grFlightPath");
  grFlightPath->SetMarkerColor(kRed);
  grFlightPath->SetMarkerStyle(6);
  // for(int i=0; i <  BaseList::getNumAbstractBases(); i++){auto gr = new TGraphAntarctica(BaseList::getAbstractBase(i), 600); gr->SetLineColor((i%10)+1);  gr->SetMarkerColor((i%10)+1); cout << gr->GetN() << endl; gr->Draw("lp");

  for(Long64_t entry=0; entry < N; entry++){

    sumtree->GetEntry(entry);
    UsefulAdu5Pat usefulPat(pat);
    for(int i=0; i <  BaseList::getNumAbstractBases(); i++){
      AntarcticCoord c = BaseList::getAbstractBase(i).getPosition(sum->realTime).as(AntarcticCoord::WGS84);
      // std::cout<< c.x<< " "<<c.y<<" "<<c.z<< std::endl;
  // Double_t getDistanceFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt); ///< Gets distance from any source in meters
      if(c.x < -60){
	       Double_t distKm = 1e-3*usefulPat.getDistanceFromSource(c.x, c.y, c.z);
      	if(distKm < maxDistKm){
      	  if(!grs[i]){
      	    grs[i] = new TGraphAntarctica();
      	    grs[i]->SetName(BaseList::getAbstractBase(i).getName());
      	    grs[i]->SetTitle(BaseList::getAbstractBase(i).getName());
      	    grs[i]->SetMarkerColor((i%7)+3);
            grs[i]->SetMarkerStyle(7);
      	  }
      	  grs[i]->SetPoint(grs[i]->GetN(), BaseList::getAbstractBase(i).getPosition(sum->realTime));
          // std::cout<<sum->realTime<<std::endl;
      	}
      }
    }

    // if((entry%1000)==0){
    grFlightPath->SetPoint(grFlightPath->GetN(), pat->longitude, pat->latitude);
    // }
  }


  for(int i=0; i < grs.size(); i++){
    if(grs[i]){
      // std::cout<< "path="<< i << std::endl;

      grs[i]->Draw("p");
    }
  }
  grFlightPath->Draw("p");

  TFile * probabilityMap = new TFile("/Users/sylarcp/anitaNeutrino/anitaBuildTool/components/UCorrelator/macro/source_maps/50_367.root");; 
  UCorrelator::ProbabilityMap * map_weighted = (UCorrelator::ProbabilityMap*) probabilityMap->Get("map_weighted");
  UCorrelator::ProbabilityMap * map_unweighted = (UCorrelator::ProbabilityMap*) probabilityMap->Get("map_unweighted");
  map_unweighted->segmentationScheme()->Draw("same",map_unweighted->getProbSums(true));
}

