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
	/** Filter that changes the response of antenna BH13 to match the rest of the antennas.  This needs to be done because any correlations with BH13 will jittery otherwise */
  class BH13Filter : public FilterOperation
  {
    public: 

      BH13Filter(){;} 

      virtual const char * tag() const { return "BH13Filter"; } 
      virtual const char * description() const { return "BH13Filter"; } 

      virtual void process(FilteredAnitaEvent * ev); 
      virtual void processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol); 
      
  }; 
	/** Filter that changes the response of antenna BH13 to back to normal, undoing the effects of BH13Filter.  This needs to be done for combining deconvolved waveforms upon which the BH13Filter has been applied */
  class AntiBH13Filter : public FilterOperation
  {
    public: 

      AntiBH13Filter(){;} 

      virtual const char * tag() const { return "AntiBH13Filter"; } 
      virtual const char * description() const { return "AntiBH13Filter"; } 

      virtual void process(FilteredAnitaEvent * ev); 
      virtual void processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol); 
      
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
