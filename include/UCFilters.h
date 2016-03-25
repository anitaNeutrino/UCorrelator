#ifndef _UCORRELATOR_FILTERS_H_
#define _UCORRELATOR_FILTERS_H_

/** This file implements all of the filters used in MyCorrelator 
 * in the new framework. 
 *
 * The AdaptiveFilter requires baselines of runs, which are handled in Baseline.h
 * */ 

#include "FilteredAnitaEvent.h" 
#include "FilterOperation.h" 
#include "TString.h" 
#include "Baseline.h"

namespace UCorrelator 
{

  /*** Apply the series of filters originally implemented in MyCorrelator to the filter strategy */ 
  void applyAbbysFilterStrategy(FilterStrategy * strategy); 

  //condition used for satellite filter
  bool antennaIsNorthFacing(FilteredAnitaEvent *ev, int trace); 


  class SimplePassBandFilter : public UniformFilterOperation
  {
    public: 
      const char * tag() const { return "SimplePassBandFilter"; }
      const char * description() const { return descStr.Data(); }

      SimplePassBandFilter(double low = 200, double high = 1200)  
        : low(low), high(high)
      {
        descStr = TString::Format("SimplePassbandFilter(%g,%g)", low,high); 
      }
    protected: 
      virtual void processOne(AnalysisWaveform *) ;

    private: 
      TString descStr; 
      double low; 
      double high;

  }; 

  class SimpleNotchFilter : public UniformFilterOperation
  {
    public: 
      SimpleNotchFilter(double minfreq, double maxfreq) 
        : min(minfreq), max(maxfreq) 
      {
        desc.Form("SimpleNotchFilter(%g,%g)",min,max); 
      }

      const char * tag() const { return "SimpleNotchFilter"; } 
      const char * description() const { return desc.Data(); } 
      virtual void processOne(AnalysisWaveform *); 
    private: 
      TString desc;  
      double min,max; 
  }; 


  class ComplicatedNotchFilter : public UniformFilterOperation 
  {
    public: 
      ComplicatedNotchFilter(double minfreqGHz, double maxfreqGHz, double temperature = 340, double gain = 75) 
        : min(minfreqGHz), max(maxfreqGHz), temperature(temperature), gain(gain) 
      {
        desc.Form("Complicated Notch Filter , T = %f K, G = %f dB", temperature, gain); 
      }

      const char * tag() const { return "ComplicatedNotchFilter"; } 
      const char * description() const { return desc.Data(); } 
      virtual  void processOne(AnalysisWaveform * aw); 
    private: 
      TString desc; ; 
      double min,max; 
      double temperature; 
      double gain; 
  }; 

  class AdaptiveFilter : public FilterOperation
  {

    public: 
      AdaptiveFilter(double dbCut, const char * baseline_dir, int navg = 1000, 
        double fmin = 0.2, double fmax = 1.2, double bandwidth = 0.026, int nFreqs = 5, double temperature = 340, double gain = 75 ) 
        : dbCut(dbCut), run(-1), navg(navg), baseline_dir(baseline_dir), baseline(0), fmin(fmin), fmax(fmax), hpol_avg(0), vpol_avg(0), 
        temperature(temperature), gain(gain)
      {
        desc_string.Form("Adaptive filter with BW %f using Baselines with %d averages and a db cut of %f\n",bw,  navg, dbCut); 
        freqs[0] = new double[nfreq]; 
        freqs[1] = new double[nfreq]; 
      }


      virtual ~AdaptiveFilter() { if (baseline) delete baseline; delete [] freqs[0]; delete [] freqs[1]; } 

      const char * tag() const { return "AdaptiveFilter"; } 
      const char * description() const { return desc_string.Data(); }
      virtual void process(FilteredAnitaEvent * event); 
      unsigned nOutputs() const { return 4 + 2*nfreq; } 
      void fillOutputs(double * vars) const; 
      const char * outputName(unsigned i) const; 

    private: 


      TString desc_string; 
      double dbCut; 
      int run; 
      int navg; 
      const char * baseline_dir; 
      Baseline * baseline; 
      double fmin, fmax; 
      double bw; 
      int nfreq;

      TGraph * hpol_avg; 
      TGraph * vpol_avg; 


      double temperature; 
      double gain; 

      //stored by print outputs
      double mean_freq[2]; 
      double strongest_cw[2]; 
      double * freqs[2]; 

  }; 



  
}



#endif
