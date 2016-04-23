#ifndef UCORRELATOR_UTIL_H
#define UCORRELATOR_UTIL_H

class UsefulAdu5Pat; 
class RawAnitaHeader; 
#include "AnitaConventions.h" 

namespace UCorrelator
{

  class AnalysisConfig; 

  double getWAISDt(const UsefulAdu5Pat * pat, const RawAnitaHeader *hdr, AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal, const AnalysisConfig * cfg = 0, double * distance = 0); 
  double getLDBDt(const UsefulAdu5Pat * pat, const RawAnitaHeader *hdr, AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal, const AnalysisConfig * cfg = 0, double * distance = 0); 

  bool isWAISVPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 
  bool isWAISHPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 

  bool isLDBVPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 
  bool isLDBHPol(const UsefulAdu5Pat * pat, const RawAnitaHeader * hdr, const AnalysisConfig * cfg = 0) ; 

}


#endif 
