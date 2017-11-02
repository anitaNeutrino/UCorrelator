#ifndef UCORRELATOR_ANALYSIS_CONFIG_H
#define UCORRELATOR_ANALYSIS_CONFIG_H

#include "AnitaConventions.h" 
class TH2; 


namespace AnitaResponse
{
class DeconvolutionMethod; 
}

namespace UCorrelator
{


  /** This class (really, it should be a struct but for CINT) stores configuration parameters that influence the behavior of th Analyzer. The defaults are defined
   * within AnalysisConfig.cc */ 
  class AnalysisConfig
  {
    public: 


      /** construct from file */ 
      AnalysisConfig(const char * config_file = 0); 
      ~AnalysisConfig();

      /** load options from file */ 
      void loadFromFile(const char * config_file); 

      unsigned correlator_nphi; /// Number of phi bins in rough correlation map 
      unsigned correlator_ntheta;  /// Number of theta bins in rough correlation map 
      double correlator_theta_lowest;  /// Lowest elevation considred, measured as positive below horizon. (negative would be above horizon) 
      double correlator_theta_highest;  /// Highest elevation considered, measured as positive above horizon. (negative would be below horizon)
      bool enable_group_delay;  /// enable group delay in interferometer 

      int zoomed_nphi; /// number of phi bins in zoomed correlation map
      int zoomed_ntheta; /// number of theta bins in zoomed correlation map
      double zoomed_dphi;  /// size of phi bins in zoomed correlation map 
      double zoomed_dtheta;  /// size of theta bins in zoomed correlation map 
      int zoomed_nant;  /// If non-zero, limit antennas considered in zoomed correlation map to nearest zoomed_nant antennas

      unsigned combine_nantennas;  /// number of antennas used to create coherent / deconvolved waveforms
      unsigned combine_npad;  /// supersampling factor for combining waveforms (i.e. how many times to pad in fourier domain. npad = 1 is super sample by 100%) 


      double saturation_threshold; /// threshold to consider a waveform saturated 

      AnitaPol::AnitaPol_t start_pol;  /// Start polarization for Analyer (kHORIZONTAL if you want just hpol or both, kVERTICAL if you want just vpol)
      AnitaPol::AnitaPol_t end_pol;  /// End polarization for Analyer (kHORIZONTAL if you want just hpo. kVERTICAL if you want just vpol or both.)

      double peak_isolation_requirement; /// Minimum distance 
      double max_peak_trigger_angle; /// Maximum distance from trigger angle to consider a peak (<=0 is any, default 0); 
      int min_peak_distance_from_unmasked; /// Minimum distance of a phi sector peak lies in from an unmasked sector. If 0, only unmasked sectors will be considered. 1 means immediate neighbors are ok. negative (default) means no check is performed.

      /** Choice of interferometer fine peak finding */ 
      enum FinePeakFindingOption_t
      {
        FinePeakFindingAbby,  ///Abby's interpolation method 
        FinePeakFindingBicubic, ///Bicubic interpolation (not implemented yet) 
        FinePeakFindingGaussianFit, ///Bivariate gaussian fit (slow) 
        FinePeakFindingQuadraticFit9, ///quadratic fit near peak, using 9 bins
        FinePeakFindingQuadraticFit16, ///quadratic fit near peak, using 16 bins
        FinePeakFindingQuadraticFit25, ///quadratic fit near peak, using 25 bins
        FinePeakFindingHistogram, ///Uses histogram means / max / rms's/correlation. Dumb but fast and foolproof. 
      } fine_peak_finding_option; 

      static const char * getPeakFindingString(FinePeakFindingOption_t opt); 

      enum ResponseOption_t
      {
        ResponseNone ,  
        ResponseSingleBRotter, /// Ben's unified respone 
        ResponseIndividualBRotter, ///Ben's individual responses 
        ResponseHarmSignalOnly, ///Harm's signal chain only thing (currently used in icemc) 
				ResponseTUFF /// response with TUFFs convolved in
      } response_option;  

      static const char * getResponseString(ResponseOption_t opt); 
      int response_npad; //number of times to pad (in freq domain) the response 


      int nmaxima; ///number of maxima computed
      bool use_bin_center; ///True to use bin center in interferometric map

      double bw_ndb; /// the bandwidth of a waveform is defined as the portion of the power spectrum  near the highest value above -ndb 
      double spectral_fit_start; 
      double spectral_fit_stop; 

      double noise_estimate_t0;  ///this is used to pick parts of the waveform for calcuating the N in SNR
      double noise_estimate_t1;  ///this is used to pick parts of the waveform for calcuating the N in SNR

      bool scale_by_cos_theta;  // Scale peak values by cos theta when picking max (due to different bin sizes) 
      bool use_offline_mask; // use offline phi masking / l1 triggers (default true)  
      double baseline_weight; //exponent for baseline weighting of correlation maps (since longer baselines give you better resolution) 


      //payload blast cuts (if 0, then ignored) 
      double max_mean_power_filtered;  
      double max_median_power_filtered; 
      double max_bottom_to_top_ratio; 


      /** Used to define when a waveform is likely from a cal pulser
       * This is only used for WAIS and Siple right now since hte LDB pulsing was so variable. 
       * */
      class Pulser
      {
        public: 
          Pulser(double offset, double distance, double dt) : GPS_offset(offset), max_distance(distance), max_dt(dt)  { ; } 
          double GPS_offset; 
          double max_distance;
          double max_dt; 
      }; 

      Pulser wais_hpol, wais_vpol, siple_hpol, siple_vpol; 

      /** For LDB pulser, we rely on a histogram generated by e.g. macros/makeLDBSelection.C for now" 
       * */
      TH2* ldb_hist() const; 
      int ldb_max_run; 

      /** When calculating stokes parameters, do you want to window the waveform around the hilbert envelope?
	  windowStokes == false: no windowing at all
	  if true:
	             stokesWindowLength =< 0: dynamic window length (default)
		     stokesWindowLength > 0: set window length
       * */
      bool windowStokes;
      int stokesWindowLength;

      /** set_bottom_first:  Re-arrange the coherently summed waveforms so that the bottom-most ring is always the
	  first antenna to be summed.
	  * */
      bool set_bottom_first;
      
      /** delay_to_center:
	  Makes it so the delays in the coherent sum are referenced to the center point */
      bool delay_to_center;

      /** TODO: this has to be loaded from file somehow */ 
      AnitaResponse::DeconvolutionMethod * deconvolution_method; 
      

      double combine_t0;
      double combine_t1;
    private: 
      mutable TH2 * the_ldb_hist; 
  };
}


#endif 
