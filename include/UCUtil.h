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


  /** Returns true if we think it's a WAIS vpol cal pulse */ 
  bool isWAISVPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 


  /** Returns true if we think it's a WAIS hpol cal pulse */ 
  bool isWAISHPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 

  /** Returns true if we think it's a LDB cal pulse */ 
  bool isLDB(const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 



  // /** Returns 0 if misses continent, 1 if it hits, 2 if it hits after adjusting theta */ 
  // int traceBackToContinent(const UsefulAdu5Pat * pat, double phi, double theta, 
  //                          double * lat, double * lon, double *alt, double * theta_adjustment_required, 
  //                          double max_theta_adjustment = 1, int max_iter = 10); 

}


#endif 
