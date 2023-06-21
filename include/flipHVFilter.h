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
  /**
   *  The purspose of this filter is to flip H and Vpol channels.
   */

  class flipHVFilter : public FilterOperation
  {
    public: 

      flipHVFilter(const char * index_file="A4ImpulseTUFFs/index.txt", bool anti = false); 
      ~flipHVFilter(); 

      virtual const char * tag() const { return "flipHVFilter"; } 
      virtual const char * description() const { return "flipHVFilter"; } 

      virtual void process(FilteredAnitaEvent * ev); 
      virtual void processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol); 

    private:
      bool anti;
      std::vector<unsigned> end_times; 
      std::vector<unsigned> indices; 
      std::vector<TGraph*> gPhase;
      std::vector<TGraph*> gMag;
      
  }; 
	/** Filter that changes the response of antenna BV13 to back to normal, undoing the effects of flipHVFilter.  This needs to be done for combining deconvolved waveforms upon which the flipHVFilter has been applied 
   * THIS IS NOW IMPLEMENTED IN TERMS OF BV13 FILTER 
   * */
  class AntiflipHVFilter : public FilterOperation
  {
    public: 

      AntiflipHVFilter(const char * index_file="A4ImpulseTUFFs/index.txt") 
        : filt(index_file,true) 
        {

        }
      ~AntiflipHVFilter() {; }

      virtual const char * tag() const { return "AntiflipHVFilter"; } 
      virtual const char * description() const { return "AntiflipHVFilter"; } 

      virtual void process(FilteredAnitaEvent * ev) { filt.process(ev); }
      virtual void processOne(AnalysisWaveform * awf, const RawAnitaHeader * header, int whichAnt, int whichPol) 
      {
        filt.processOne(awf,header,whichAnt,whichPol); 
      }
     
    private:
      flipHVFilter filt; 
  }; 
//  class timePadFilter : public FilterOperation
//  {
//    public: 
//
//      timePadFilter(int samples):fSamples(samples){;} 
//
//      virtual const char * tag() const { return "flipHVFilter"; } 
//      virtual const char * description() const { return "flipHVFilter"; } 
//
//      virtual void process(FilteredAnitaEvent * ev); 
//		private:
//			int fSamples;
//     
//  }; 

  /** Filter to attempt to convert A3 events to A4 events by convolving in the TUFF response for ~direct comparison between the two flights. Config argument is based on the config letter from Oindree's elog 711.  Default is the most common config. **/

//  class A3toA4ConversionFilter : public FilterOperation
//  {
//    public:
//      A3toA4ConversionFilter(char config='B');
//      ~A3toA4ConversionFilter();
//
//      virtual const char * tag() const { return "A3toA4ConversionFilter"; }
//      virtual const char * description() const { return "A3toA4ConversionFilter"; }
//
//      virtual void process(FilteredAnitaEvent * ev);
//      virtual void processOne(AnalysisWaveform * awf, const RawAnitaHeader * header=0, int whichAnt=0, int whichPol=0);
//
//    private:
//      TGraph* gReal;
//      TGraph* gImag;
//      double phaseShift;
//  };


 
}



#endif
