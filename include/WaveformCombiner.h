#ifndef UCORRELATOR_WAVEFORM_COMBINER_H
#define UCORRELATOR_WAVEFORM_COMBINER_H

/** This combines multiple waveforms. It uses supersampling and interpolation, but could in principle use a Fourier phase shift 
 *
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 **/

#include "AnalysisWaveform.h"
#include "AnitaConventions.h"
#include <vector> 


class FilteredAnitaEvent; 

namespace AnitaResponse
{
  class ResponseManager; 

}

namespace UCorrelator
{

  class WaveformCombiner
  {

    public: 

    WaveformCombiner(int nantennas = 10, int npad = 3, bool useUnfiltered = true, bool deconvolve = false, const AnitaResponse::ResponseManager * response = 0, bool enableALFAHack=true); 

      /** Sets the responses used when deconvolving. 
       * If a response is 0, nothing is done. 
       * Responses may be set globally, globally for polarizations or per antenna 
       */ 

    void setResponseManager(const AnitaResponse::ResponseManager * rm)  { responses = rm;} 

      virtual ~WaveformCombiner(); 

      /** Combines the waveforms from the given event */ 
      void combine(double phi, double theta, const FilteredAnitaEvent * event, AnitaPol::AnitaPol_t pol, ULong64_t disallowed = 0); 

      const AnalysisWaveform * getCoherent() const { return &coherent; }
      const AnalysisWaveform * getDeconvolved() const; 
      const TGraph * getCoherentAvgSpectrum() const { return &coherent_avg_spectrum; } 
      const TGraph * getDeconvolvedAvgSpectrum() const { return &deconvolved_avg_spectrum; } 
      void setNPad(int new_npad) { npad = new_npad; }
      void setNAntennas(int n) { nant = n; antennas.resize(n);  }
      void setDeconvolve(bool deconvolve) {do_deconvolution = deconvolve ;}
      void setUseUnfiltered(bool raw_opt) {use_raw = raw_opt;}
      void setGroupDelayFlag(bool opt) { enable_group_delay = opt; } 
      bool wasAlfaFiltered() { return alfa_hack; } 

      int getNAntennas() const { return nant; }  
      const int * getUsedAntennas() const { return &antennas[0]; } 

      /* setBottomFirst():
       getClosestAntennas() returns the closest antenna in phi, which could be on any ring.
    	 this will force the bottom antenna (largest ant#) to be the seed waveform so things don't jump as much
   	   default is false (like it was before)
      */
      void setBottomFirst(bool opt) { bottom_first = opt; }

      /* setDelayToCenter();
	    changes the delays in combining from being to the first antenna in the array, to being towards the centerpoint
	    of the instrument (0,0,0) */
    	void setDelayToCenter(bool opt) {delay_to_center = opt; }

      /** Static helper used to combine arbitrary waveforms */
      static AnalysisWaveform *  combineWaveforms(int nwf, const AnalysisWaveform * wfs, const double * delays, const double * scales = 0, AnalysisWaveform * output = 0); 

    private: 

      
      AnalysisWaveform coherent; 
      AnalysisWaveform deconvolved; 
      TGraphAligned coherent_avg_spectrum; 
      TGraphAligned deconvolved_avg_spectrum; 


      const AnalysisWaveform * wf(const FilteredAnitaEvent*, int ant, AnitaPol::AnitaPol_t pol); 
      const AnitaResponse::ResponseManager * responses; 
      int npad; 
      int nant; 
      bool use_raw; 
      bool do_deconvolution; 
      bool enable_group_delay; 
      bool alfa_hack; 
      bool bottom_first;
      bool delay_to_center;
      std::vector<int> antennas; 
  };

}

#endif
