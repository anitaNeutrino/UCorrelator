#ifndef _UCORRELATOR_SINE_SUB_CACHE_H 
#define _UCORRELATOR_SINE_SUB_CACHE_H

/* 
   Ben Strutt is writing a UCorrelator class?
 */
#include "SineSubtract.h"
#include "UCFilters.h"
#include "AnitaConventions.h"

class TFile;
class TTree;

namespace UCorrelator {


/** 
 * @class Handles loading/saving of sine subtract results since it's expensive
 */
class SineSubtractCache {

 public:

  static void makeCache(int run, SineSubtractFilter* ssf); // generation function
  
  SineSubtractCache(const char* descr);
  virtual ~SineSubtractCache();
  const FFTtools::SineSubtractResult* getResult(UInt_t eventNumber, AnitaPol::AnitaPol_t pol, Int_t antenna);


  Bool_t fDebug;
 private:

  static TString branchName(AnitaPol::AnitaPol_t pol, Int_t ant);
  static TString fileName(const char* specDir, UInt_t hash, Int_t run);

  void loadRun(int run);
  
  TFile* fFile;
  TTree* fTree;  
  UInt_t fDescHash; 
  TString fSpecDir;
  Int_t fCurrentRun;
  UInt_t fLastEventNumber;
  FFTtools::SineSubtractResult* results[AnitaPol::kNotAPol][NUM_SEAVEYS];

};
  


}



#endif
