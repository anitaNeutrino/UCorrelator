#ifndef _UCORRELATOR_SPECTRUM_PARAMETERS_H
#define _UCORRELATOR_SPECTRUM_PARAMETERS_H

class TGraph; 

#include "AnitaEventSummary.h" 

namespace UCorrelator
{
  class AnalysisConfig; 
  namespace spectrum
  {
    void fillSpectrumParameters( const TGraph * spectrum, const TGraph * average, AnitaEventSummary::WaveformInfo * winfo, const AnalysisConfig * config); 


  }
}



#endif
