#ifndef _BH13FILTER_H_
#define _BH13FILTER_H_


#include "FilteredAnitaEvent.h" 
#include "FilterOperation.h" 
#include "TString.h" 
#include "DigitalFilter.h" 
#include "Baseline.h"
#include "SineSubtract.h" 
#include "AnalysisWaveform.h" 

namespace UCorrelator 
{

  class BH13Filter : public FilterOperation
  {
    public: 

      BH13Filter(){;} 

      virtual const char * tag() const { return "BH13Filter"; } 
      virtual const char * description() const { return "BH13Filter"; } 

      virtual void process(FilteredAnitaEvent * ev); 
      
  }; 
  class timePadFilter : public FilterOperation
  {
    public: 

      timePadFilter(int samples):fSamples(samples){;} 

      virtual const char * tag() const { return "BH13Filter"; } 
      virtual const char * description() const { return "BH13Filter"; } 

      virtual void process(FilteredAnitaEvent * ev); 
		private:
			int fSamples;
     
  }; 

 
}



#endif
