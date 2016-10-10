#ifndef UCORRELATOR_WAVEFORM_COMBINER_H
#define UCORRELATOR_WAVEFORM_COMBINER_H

/** This combines multiple waveforms. It uses supersampling and interpolation, but could in principle use a Fourier phase shift 
 *
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 **/

#include "AnalysisWaveform.h"
#include "AnitaConventions.h"
#include <vector> 
#include <stdint.h>


class FilteredAnitaEvent; 

namespace UCorrelator
{
  class ResponseManager; 

  class WaveformCombiner
  {

    public: 

      WaveformCombiner(int nantennas = 10, int npad = 3, bool useUnfiltered = true, bool deconvolve = false, const ResponseManager * response = 0); 

      /** Sets the responses used when deconvolving. 
       * If a response is 0, nothing is done. 
       * Responses may be set globally, globally for polarizations or per antenna 
       */ 

      void setResponseManager(const ResponseManager * rm)  { responses = rm;} 

      virtual ~WaveformCombiner(); 

      /** Combines the waveforms from the given event */ 
      void combine(double phi, double theta, const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol, uint64_t disallowed = 0); 

      const AnalysisWaveform * getCoherent() const { return &coherent; }
      const AnalysisWaveform * getDeconvolved() const; 
      const TGraph * getCoherentAvgSpectrum() const { return &coherent_avg_spectrum; } 
      const TGraph * getDeconvolvedAvgSpectrum() const { return &deconvolved_avg_spectrum; } 
      void setNPad(int new_npad) { npad = new_npad; }
      void setNAntennas(int n) { nant = n; }
      void setDeconvolve(bool deconvolve) {do_deconvolution = deconvolve ;}
      void setUseUnfiltered(bool raw_opt) {use_raw = raw_opt;}
      void setGroupDelayFlag(bool opt) { enable_group_delay = opt; } 

      /** Static helper used to combine arbitrary waveforms */
      static AnalysisWaveform *  combineWaveforms(int nwf, const AnalysisWaveform * wfs, const double * delays, const double * scales = 0, AnalysisWaveform * output = 0); 

    private: 

      
      AnalysisWaveform coherent; 
      AnalysisWaveform deconvolved; 
      TGraphAligned coherent_avg_spectrum; 
      TGraphAligned deconvolved_avg_spectrum; 


      const ResponseManager * responses; 
      int npad; 
      int nant; 
      bool use_raw; 
      bool do_deconvolution; 
      bool enable_group_delay; 

  };

}

#endif
