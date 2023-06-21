#ifndef _FLIPHVFILTER_H_
#define _FLIPHVFILTER_H_


#include <vector>
#include "FilteredAnitaEvent.h" 
#include "FilterOperation.h" 
#include "TString.h" 
#include "DigitalFilter.h" 
#include "Baseline.h"
#include "SineSubtract.h" 
#include "AnalysisWaveform.h" 


namespace UCorrelator 
{
  /*
   *  The purspose of this filter is to flip H and Vpol channels.
   */

  class flipHVFilter : public FilterOperation {
  
    public: 

      flipHVFilter(); 
      ~flipHVFilter(); 

      virtual const char * tag() const { return "flipHVFilter"; } 
      virtual const char * description() const { return "flipHVFilter"; } 

      virtual void process(FilteredAnitaEvent * ev); 
      virtual void processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol);
  };  
}


#endif
