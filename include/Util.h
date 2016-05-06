#ifndef UCORRELATOR_UTIL_H
#define UCORRELATOR_UTIL_H

class UsefulAdu5Pat; 
class RawAnitaHeader; 
#include "AnitaConventions.h" 

/** \file 
 * random stuff goes here */ 

namespace UCorrelator
{

  class AnalysisConfig; 

  /** Get difference in time between trigger time and expected arrival of a WAIS cal pulse */
  double getWAISDt(const UsefulAdu5Pat * pat, const RawAnitaHeader *hdr, AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal, const AnalysisConfig * cfg = 0, double * distance = 0); 

  /** Get difference in time between trigger time and expected arrival of an LDB cal pulse */
  double getLDBDt(const UsefulAdu5Pat * pat, const RawAnitaHeader *hdr, AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal, const AnalysisConfig * cfg = 0, double * distance = 0); 


  /** Returns true if we think it's a WAIS vpol cal pulse */ 
  bool isWAISVPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 


  /** Returns true if we think it's a WAIS hpol cal pulse */ 
  bool isWAISHPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 

  /** Returns true if we think it's a LDB vpol cal pulse */ 
  bool isLDBVPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 

  /** Returns true if we think it's a LDB hpol cal pulse */ 
  bool isLDBHPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 

}


#endif 
