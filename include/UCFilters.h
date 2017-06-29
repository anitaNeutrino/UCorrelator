#ifndef _UCORRELATOR_FILTERS_H_
#define _UCORRELATOR_FILTERS_H_

/** \file This file implements all of the filters used in MyCorrelator 
 * in the new framework. 
 * */ 

#include "FilterOperation.h" 
#include "TString.h" 
#include "DigitalFilter.h" 
#include "Baseline.h"
#include "SineSubtract.h" 
#include "AnalysisWaveform.h" 

class FilteredAnitaEvent; 
class RawAnitaHeader; 


namespace AnitaResponse{
  class ResponseManager; 
  class DeconvolutionMethod; 
}

namespace UCorrelator 
{

  class SpectrumAverageLoader; 


  /** Convenience function for getting a FilterStrategy with a given string
   * key. 
   *
   *
   *  !!! all of these include an ALFA filter at the very end for A3 !!! 
   *
   *  They may be chained together as well (require a + to separate, e.g.
   *  decon+sinsub_05_0) 
   *
   * Valid keys so far are: 
   *
   *  sinsub_%02d_%d[_option]  sin subtraction format: (100 * minimum power
   *  reduction, num failed iter) 
   *
   *  Options include: 
   *
   *     _ad_%d Make adaptive (with peakiness exponent). peakiness exponent always between 1 and 10, so if e.g. is 15, interpreted as 1.5
   *
   *     Example: sinsub_05_0_ad_2 
   *       
   *     _env_%s[_%g_%g]  Use envelope detection. %s should be one of
   *     "hilbert","rms" or "peak" 
   *
   *     If you pass additional parameters they are passed as parameters to 
   *     setEnvelopeOption(otherwise the defaults are used). 
   *
   *     _peakfind_%s[_%g_%g_%g_%g]  Modify peak-finding preferences 
   *
   *     %s should be one of "global", "neighbor," "tspectrum",  or "savgolsub" 
   *
   *      If you pass additional parameters, they are passed as params to
   *      setPeakFindingOption (otherwise defaults used)
   *
   * 
   *  adsinsub_%d_%02d_%d  *  adaptive sin subtraction
   *  Shortcut for to sinsub_%02d_%d_ad_%d  (where the first %d becomes the last %d) 
   *
   *  butter_%d_%d adaptive butterworth
   *    format : (peakiness threshold,order)
   * 
   *    peakiness threshold is interpreted as necessarily between 1 and 10, so if
   *    it's e.g. 15, it becomes 1.5 
   *      
   *  minphase_%d adaptive minimum phase 
   *    format : (peakiness exponent)
   *
   *    peakiness exponent is interpreted as necessarily between 1 and 10, so if
   *    it's e.g. 15, it becomes 1.5 
   *
   *
   *  brickwall_%d_%d  brick wall with peakiness threshold
   *    format (peakiness threshold, fill notch)
   *
   *   peakiness threshold is interpreted as necessarily
   *   between 1 and 10, so if it's e.g. 15, it becomes 1.5 
   *   if fill notch not 0,  notch is filled according to spectrum average
   *         
 
   *  decon deconvolve filter  (using ``best'' known response) . Should be used in conjunction with others. 
   *
   *  geom  geometric filter ( preliminary support) 
   *
   *  abby  attempt at duplicating Abby's filter strategy (probably doesn't work for A3) 
   *
   *  Returns a human readable description.  
   */ 
  const char * fillStrategyWithKey(FilterStrategy * fillme, const char * key); 

  /* same as above, but returns a newly allocated strategy and no description */ 
  FilterStrategy * getStrategyWithKey(const char * key) ; 
  


  /** Condition used for satellite filter */ 
  bool antennaIsNorthFacing(FilteredAnitaEvent *ev, int ant, AnitaPol::AnitaPol_t); 


  /** 
   * This filter is a brick-wall notch filter that then fills in the notch using thermal noise. 
   */
  class ComplicatedNotchFilter : public UniformFilterOperation 
  {
    public: 
      /** Build a ComplicatedNotchFilter. The temperature and gain are used to model the thermal noise.  */
      ComplicatedNotchFilter(double minfreqGHz, double maxfreqGHz, double temperature = 170, double gain = 65) 
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



  /** Attempted reimplementation of Abby's Adaptive Filter. Dont' know if it works fully yet... */
  class AdaptiveFilterAbby : public FilterOperation
  {

    public: 

      AdaptiveFilterAbby(double dbCut, const char * baseline_dir, int navg = 1000, 
        double fmin = 0.2, double fmax = 1.2, double bandwidth = 0.026, int nFreqs = 5, double temperature = 340, double gain = 75 ) 
        : dbCut(dbCut), run(-1), navg(navg), baseline_dir(baseline_dir), baseline(0), fmin(fmin), fmax(fmax), bw(bandwidth), nfreq(nFreqs), hpol_avg(0), vpol_avg(0), 
        temperature(temperature), gain(gain)
      {
        desc_string.Form("Abby's adaptive filter with BW %f using Baselines with %d averages and a db cut of %f\n",bw,  navg, dbCut); 
        freqs[0] = new double[nfreq]; 
        freqs[1] = new double[nfreq]; 
      }


      virtual ~AdaptiveFilterAbby() { if (baseline) delete baseline; delete [] freqs[0]; delete [] freqs[1]; } 

      const char * tag() const { return "AdaptiveFilter"; } 
      const char * description() const { return desc_string.Data(); }
      virtual void process(FilteredAnitaEvent * event); 
      unsigned nOutputs() const { return 6 ; } 
      virtual void fillOutput(unsigned i, double * vars) const; 
      const char * outputName(unsigned i) const; 
      unsigned outputLength(unsigned i) const { return i < 4 ? 1 : nfreq; } 
      const TGraph * getHpolAvg() const { return hpol_avg; }
      const TGraph * getVpolAvg() const { return vpol_avg; }
      const TGraph * getLastPowerHpol() const { return &powers[0]; }
      const TGraph * getLastPowerVpol() const { return &powers[1]; }

    private: 


      TString desc_string; 
      double dbCut; 
      TGraph powers[2]; 
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




  class SineSubtractFilter
    : public FilterOperation
  {
    public: 
      SineSubtractFilter(double min_power_ratio = 0.05, int max_failed_iter = 0,  int nfreq_bands = 0, const double *  freq_bands_start = 0, const double * freq_bands_end = 0, int nstored_freqs = 5); 

      /** Make the filter adaptive using a SpectrumAverageLoader. If null passed, adaptiveness turned off.  */ 
      void makeAdaptive(const SpectrumAverageLoader *avg = 0, double peakiness_exp = 1); 

      virtual ~SineSubtractFilter();  
      void setInteractive(bool set); 

      const char * tag() const { return "SineSubtractFilter"; }
      const char * description() const { return desc_string.Data(); } 
      unsigned nOutputs() const { return 1+2*(2+3*nstored_freqs); }
      const char * outputName(unsigned i) const{ return output_names[i].Data(); }
      unsigned outputLength(unsigned i) const{ return i == 0 ? 1 : NUM_SEAVEYS; } 
      void fillOutput(unsigned i, double * vars) const; 
      void setVerbose(bool v); 

      void setPeakFindingOption(FFTtools::SineSubtract::PeakFindingOption option, const double * params = 0);
      void setEnvelopeOption(FFTtools::SineSubtract::EnvelopeOption option, const double * params = 0);
      virtual void process(FilteredAnitaEvent * ev); 
      const FFTtools::SineSubtract* sinsub(AnitaPol::AnitaPol_t pol, int ant) const { return subs[pol][ant] ;} 
    private:
      void processOne(AnalysisWaveform* wf,const RawAnitaHeader* h,int i, int pol); 
      FFTtools::SineSubtract * subs[2][NUM_SEAVEYS]; 
      double min_power_ratio; 
      const SpectrumAverageLoader * spec; 
      TGraph * reduction[2][NUM_SEAVEYS]; 
      unsigned last_t; 
      TString desc_string; 
      std::vector<TString> output_names;
      int nstored_freqs; 
      double adaptive_exp; 
      int max_failed; 
  };

  class AdaptiveBrickWallFilter : public FilterOperation
  {
    public: 
      AdaptiveBrickWallFilter(const UCorrelator::SpectrumAverageLoader * spec, double thresh=2, bool fillNotch = true);  

      const char * tag() const { return "AdaptiveBrickWallFilter"; } 
      const char * description() const{ return desc_string.Data(); } 
      virtual void process(FilteredAnitaEvent *ev); 
      virtual ~AdaptiveBrickWallFilter();
    private:
      TString desc_string; 
      const SpectrumAverageLoader * avg; 
      double threshold; 
      bool fill; 
      int last_bin; 
      int instance; 
      TH1* sp[2][NUM_SEAVEYS]; 
  }; 



  class AdaptiveMinimumPhaseFilter : public FilterOperation
  {

    public: 
      AdaptiveMinimumPhaseFilter(const SpectrumAverageLoader * avg, double exponent = -2, int npad =3); 

      const char * tag() const { return "AdaptiveMinimumPhaseFilter"; } 
      const char * description() const{ return desc_string.Data(); } 
      virtual void process(FilteredAnitaEvent *ev); 
      virtual ~AdaptiveMinimumPhaseFilter();  
      TGraph * getCurrentFilterTimeDomain(AnitaPol::AnitaPol_t pol, int i) const; 
      TGraph * getCurrentFilterPower(AnitaPol::AnitaPol_t pol, int i) const; 

    private: 
      TString desc_string; 
      const SpectrumAverageLoader * avg; 
      int npad; 
      double exponent; 
      int last_bin; 
      FFTWComplex * filt[2][NUM_SEAVEYS]; 
      int size[2][NUM_SEAVEYS]; 
  }; 


  class AdaptiveButterworthFilter : public FilterOperation 
  {
    public: 
      AdaptiveButterworthFilter(const SpectrumAverageLoader *avg, double peakiness_threshold = 2, int order = 2, double width = 0.05) ; 

      virtual void process(FilteredAnitaEvent *ev) ; 
      virtual ~AdaptiveButterworthFilter() {; } 
      const char * tag() const { return "AdaptiveButterworthFilter"; } 
      const char * description() const{ return desc_string.Data(); } 
      const FFTtools::DigitalFilterSeries * getFilter(AnitaPol::AnitaPol_t pol, int ant) const { return &filters[pol][ant]; } 

    private: 
      TString desc_string; 
      const SpectrumAverageLoader * avg; 
      double threshold; 
      int last_bin; 
      int order; 
      double width; 
      FFTtools::DigitalFilterSeries filters[2][NUM_SEAVEYS]; 
  }; 



  class CombinedSineSubtractFilter
    : public FilterOperation
  {
    public: 
      CombinedSineSubtractFilter(double min_power_ratio = 0.05, int max_failed_iter = 0, int nfreq_bands = 0, const double *  freq_bands_start = 0, const double * freq_bands_end = 0, int nstored_freqs = 5); 

      virtual ~CombinedSineSubtractFilter();  

      const char * tag() const { return "CombinedSineSubtractFilter"; }
      const char * description() const { return desc_string.Data(); } 
      unsigned nOutputs() const { return 1+(2+3*nstored_freqs); }
      const char * outputName(unsigned i) const{ return output_names[i].Data(); }
      unsigned outputLength(unsigned i) const; 
      void fillOutput(unsigned i, double * vars) const; 
      virtual void process(FilteredAnitaEvent * ev); 
      const FFTtools::SineSubtract* sinsub(AnitaPol::AnitaPol_t pol, int phi) const { return subs[phi][pol] ;} 
      void setInteractive(bool set); 
    private:
      FFTtools::SineSubtract * subs[NUM_PHI][2]; 
      TString desc_string; 
      std::vector<TString> output_names;
      int nstored_freqs; 
  };
}



#endif
