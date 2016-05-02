#ifndef UCORRELATOR_WAVEFORM_COMBINER_H
#define UCORRELATOR_WAVEFORM_COMBINER_H


#include "AnalysisWaveform.h"
#include "AnitaConventions.h"
class FilteredAnitaEvent; 
#include <stdint.h>


namespace UCorrelator
{
  class AbstractResponse; 

  class WaveformCombiner
  {

    public: 
      WaveformCombiner(int nantennas = 10, int npad = 3, bool useUnfiltered = true, bool deconvolve = false, const char * responseDir = ""); 

      /** Sets the responses used when deconvolving. 
       * If a response is 0, nothing is done. 
       * Responses may be set globally, globally for polarizations or per antenna 
       */ 
      void setResponse(const AbstractResponse * response = 0, AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol, int antenna = -1); 
      void loadResponsesFromDir(const char * dir); 

      virtual ~WaveformCombiner(); 

      /** Combines the waveforms from the given event */ 
      void combine(double phi, double theta, const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol, uint64_t disallowed = 0); 
      const AnalysisWaveform * getCoherent() const { return &coherent; }
      const AnalysisWaveform * getDeconvolved() const; 
      void setNPad(int new_npad) { npad = new_npad; }
      void setNAntennas(int n) { nant = n; }
      void setDeconvolve(bool deconvolve) {do_deconvolution = deconvolve ;}
      void setUseUnfiltered(bool raw_opt) {use_raw = raw_opt;}
      void setGroupDelayFlag(bool opt) { enable_group_delay = opt; } 

      static AnalysisWaveform *  combineWaveforms(int nwf, const AnalysisWaveform * wfs, const double * delays, const double * scales = 0, AnalysisWaveform * output = 0); 

    private: 

      
      AnalysisWaveform coherent; 
      AnalysisWaveform deconvolved; 


      const AbstractResponse * responses[2][NUM_SEAVEYS]; 
      int npad; 
      int nant; 
      bool use_raw; 
      bool do_deconvolution; 
      bool enable_group_delay; 

  };

}

#endif
