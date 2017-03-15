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
#include "BasicFilters.h" 
#include "AnalysisConfig.h"
#include "AnitaDataset.h"
#include "RawAnitaHeader.h"




int min_ldb = 130; 
int max_ldb = 164; 

int min_wais = 332; 
int max_wais = 354; 


std::vector<FilterStrategy*> strategies; 
std::vector<std::string> names; 

void addStrategy( FilterStrategy * s, const char * name) { strategies.push_back(s); names.push_back(name); } 



/** add filter strategies here ! */ 

void setupFilters(TFile* out) 
{

  /** Sine subtraction */ 
  FilterStrategy * sinsub= new FilterStrategy; 
  double fmins[1] = {0.18}; 
  double fmaxs[1] = {1.3}; 
  sinsub->addOperation(new UCorrelator::SineSubtractFilter(0.05, 0, 4,1,fmins,fmaxs)); 
  sinsub->addOperation(new ALFAFilter); 
  addStrategy(sinsub, "sinsub"); 


  /** Peakiness Detecting Sine Subtraction */ 



  /** Adaptive Butterworth Filter */ 



  /** Minimum Phase **/ 

}




int main(int nargs, char ** args) 
{

  FFTtools::loadWisdom("wisdom.dat"); 
  int run = atoi(args[1]); 

  bool isWAIS = run <=min_wais && run >= max_wais; 
  bool isLDB = run <=min_ldb && run >= max_ldb; 
  bool isBG  = !isWAIS && !isLDB; 

  AnitaDataset d(run, isBG); //only use decimated if background 

  int max = nargs > 2 ? atoi(args[2]) : 0; 

  TString outname; 
  const char * label = isWAIS ? "wais" : isLDB ? "ldb" : "bg"; 

  if (max) outname.Form("filter/%d_%s_max_%d.root", run, label, max); 
  else outname.Form("filter/%d_%s_max_%d.root", run, label, max); 

  TFile ofile(outname, "RECREATE"); 
  setupFilters(&ofile); 


  if (!ofile.IsOpen())
  {
    return 1; 
  }

  UCorrelator::AnalysisConfig cfg; 
  UCorrelator::Analyzer analyzer(&cfg); 

  AnitaEventSummary * sum = new AnitaEventSummary; 
  RawAnitaHeader *hdr = 0 ; 
  UsefulAdu5Pat *patptr = 0; 

  std::vector<TTree*> trees; 

  TTree * friendly = new TTree("aux", "Auxdata (headers/gps)"); 
  friendly->Branch("header",&hdr); 
  friendly->Branch("pat",&patptr); 

  for (size_t i = 0; i < strategies.size(); i++) 
  {
    TString tname; tname.Form("filter %s (on %s)",names[i].c_str(),label);
    trees[i] = new TTree(names[i].c_str(),tname.Data()); 
    trees[i]->Branch("summary",sum); 
    trees[i]->AddFriend(friendly); 
  }


  int ndone = 0; 

  for (int i = 0; i < d.N(); i++) 
  {
    d.getEntry(i); 

    hdr = d.header(); 
    if (isLDB) 
    {
      if (UCorrelator::isLDB(d.header()))
      {
         printf("----(LDB event %d ( idx=%d))-----\n",hdr->eventNumber,i); 
      }
      else continue; 
    }

    UsefulAdu5Pat pat(d.gps()); 

    if (isWAIS) 
    {
      if (UCorrelator::isWAISHPol(&pat, d.header()))
      {
         printf("----(WAIS event %d (idx=%d))-----\n",hdr->eventNumber,i); 
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

    //loop over strategies 
    for (size_t s = 0; s < strategies.size(); s++) 
    {
      FilteredAnitaEvent ev(d.useful(), strategies[s], d.gps(), d.header()); 
      analyzer.analyze(&ev, sum); 
      ofile.cd(); 
      trees[s]->Fill(); 
    }

    if (max && ndone > max) break; 
    ndone++; 
  }


  ofile.cd(); 

//  friendly->Write(); 
//  for (size_t i = 0; i < trees.size(); i++) trees[i]-Write(); 

  FFTtools::saveWisdom("wisdom.dat"); 
  ofile.Write(); 
}




