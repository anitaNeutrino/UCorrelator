#ifndef UCORRELATOR_ANALYSIS_CONFIG_H
#define UCORRELATOR_ANALYSIS_CONFIG_H

#include "AnitaConventions.h" 

namespace UCorrelator
{

  class AnalysisConfig
  {
    public: 
      AnalysisConfig(); 
#ifdef ENABLE_LIBCONFIG
      AnalysisConfig(const char * config_file = 0); 
      void loadFromFile(const char * config_file); 
#endif


      unsigned correlator_nphi; 
      unsigned correlator_ntheta; 
      double correlator_theta_lowest; 
      double correlator_theta_highest; 
      bool enable_group_delay; 

      int zoomed_nphi; 
      int zoomed_ntheta; 
      double zoomed_dphi; 
      double zoomed_dtheta; 
      int zoomed_nant; 

      unsigned combine_nantennas; 
      unsigned combine_npad; 


      double saturation_threshold; 

      AnitaPol::AnitaPol_t start_pol; 
      AnitaPol::AnitaPol_t end_pol; 

      double peak_isolation_requirement; 

      enum FinePeakFindingOption_t
      {
        FinePeakFindingAbby,  //Abby's interpolation method 
        FinePeakFindingBicubic, //Bicubic interpolation
        FinePeakFindingGaussianFit, //Bivariate gaussian fit
        FinePeakFindingQuadraticFit9, //quadratic fit near peak, using 9 bins
        FinePeakFindingQuadraticFit16, //quadratic fit near peak, using 16 bins
        FinePeakFindingQuadraticFit25, //quadratic fit near peak, using 25 bins
      } fine_peak_finding_option; 


      int nmaxima; //number of maxima computed
      bool use_bin_center; 

      double bw_ndb; 
      double noise_estimate_t0; 
      double noise_estimate_t1; 


      class Pulser
      {
        public: 
          Pulser(double offset, double distance, double dt) : GPS_offset(offset), max_distance(distance), max_dt(dt)  { ; } 
          double GPS_offset; 
          double max_distance;
          double max_dt; 
      }; 

      Pulser wais_hpol, wais_vpol, ldb_hpol, ldb_vpol;  
  };
}


#endif 
